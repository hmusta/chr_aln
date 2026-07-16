#pragma once

#include "helpers.hpp"
#include "wfa_main.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>

#include <cstdint>
#include <cstddef>

#include "monotone_diagonal_profile.hpp"

template <bool inv_open_temp = false, bool inv_ext_temp = false, bool is_fwd_temp = true>
class WFAIterator {
  public:
    static constexpr bool inv_open = inv_open_temp;
    static constexpr bool inv_ext = inv_ext_temp;
    static constexpr bool is_fwd = is_fwd_temp;

    using Container = chr_aln::ForwardDiagonalProfile;

    using const_iterator = typename Container::const_iterator;
    using string_const_iterator = std::conditional_t<
        is_fwd,
        std::string_view::const_iterator,
        std::reverse_iterator<std::string_view::const_iterator>
    >;

    using Callback = std::function<void(Offset, Offset, Offset)>;

    static Diag get_diag(SOffset q, SOffset t) { return t - q; }

    WFAIterator(const ScoreModel &score_model,
                std::string_view::const_iterator qbegin,
                std::string_view::const_iterator qend,
                std::string_view::const_iterator tbegin,
                std::string_view::const_iterator tend,
                UDiag min_wavefront_length = max_diag_width,
                SOffset max_distance_threshold = max_offset)
          : score_model_(score_model),
            qsize_(std::distance(qbegin, qend)),
            tsize_(std::distance(tbegin, tend)),
            empty_counter_(0),
            p_(0),
            min_k_(get_diag(-1, tsize_ + 1)),
            max_k_(get_diag(qsize_ + 1, -1)),
            min_q_(qsize_ + 1),
            max_q_(-1),
            min_t_(tsize_ + 1),
            max_t_(-1),
            frp_q_(-1),
            frp_t_(-1),
            frp_p_(0),
            frp_ext_(0),
            endpoints_(tsize_ + 2 + qsize_ + 2),
            qpoints_(tsize_ + 2 + qsize_ + 2, tsize_ + 2 + qsize_ + 2),
            tpoints_(tsize_ + 2 + qsize_ + 2, tsize_ + 2 + qsize_ + 2),
            global_min_(endpoints_.size()),
            min_wf_len_(min_wavefront_length),
            max_dist_(max_distance_threshold) {
        if constexpr (is_fwd) {
            qbegin_ = qbegin;
            qend_ = qend;
            tbegin_ = tbegin;
            tend_ = tend;
        } else {
            qbegin_ = std::make_reverse_iterator(qend);
            qend_ = std::make_reverse_iterator(qbegin);
            tbegin_ = std::make_reverse_iterator(tend);
            tend_ = std::make_reverse_iterator(tbegin);
        }
        assert(std::distance(qbegin_, qend_) == static_cast<ssize_t>(qsize_));
        assert(std::distance(tbegin_, tend_) == static_cast<ssize_t>(tsize_));

        std::get<0>(table_.emplace_back().emplace_back(0, def_elem).second) = 0;
        process();

        assert(min_k_ == 0);
        assert(max_k_ == 0);
    }

    std::string_view get_target() const {
        if constexpr (is_fwd) {
            return std::string_view(tbegin_, std::distance(tbegin_, tend_));
        } else {
            return std::string_view(tend_.base(), std::distance(tbegin_, tend_));
        }
    }

    std::string_view get_query() const {
        if constexpr (is_fwd) {
            return std::string_view(qbegin_, std::distance(qbegin_, qend_));
        } else {
            return std::string_view(qend_.base(), std::distance(qbegin_, qend_));
        }
    }

    SOffset max_q() const { return qsize_; }
    SOffset max_t() const { return tsize_; }

    Diag get_min_diag() const { return min_k_; }
    Diag get_max_diag() const { return max_k_; }

    SOffset get_global_begin() const { return min_q_; }
    SOffset get_global_end() const { return max_q_; }
    SOffset get_global_r_begin() const { return min_t_; }
    SOffset get_global_r_end() const { return max_t_; }
    Penalty get_p() const { return p_; }

    bool empty() const { return empty_counter_ >= table_.alloc_size(); }

    std::tuple<SOffset, SOffset, Penalty, SOffset> get_frp() const {
        return std::make_tuple(frp_q_, frp_t_, frp_p_, frp_ext_);
    }

    UDiag get_max_antidiag() const { return frp_q_ + frp_t_; }
    bool end_reached() const { return get_max_antidiag() == tsize_ + qsize_; }

    const Container& diag_data(Diag diag) const { return endpoints_[diag_to_idx(diag)]; }
    const_iterator begin(Diag diag) const { return diag_data(diag).begin(); }
    const_iterator end(Diag diag) const { return diag_data(diag).end(); }

    std::pair<const_iterator, const_iterator> get_min(Diag diag, Offset q_b, Offset q_e) const {
        return diag_data(diag).query_min(q_b, q_e);
    }

    std::pair<const_iterator, const_iterator> get_min_range(Diag diag, Offset q_b, Offset q_e) const {
        return diag_data(diag).query_min_range(q_b, q_e);
    }

    template <bool observed_min = false, bool observed_max = false>
    std::pair<Diag, Diag> get_diag_range_from_q(SOffset q_b, SOffset q_e) {
        // k = t - q => t = k + q
        // we need t >= t_obs_min and t <= t_obs_max
        // k + q >= t_obs_min => k >= t_obs_min - q
        // k + q <= tsize_ => k <= tsize_ - q
        Diag d_front = std::max<Diag>(get_min_diag(), (observed_min ? min_t_ : 0) - q_b);
        assert(d_front + q_b >= (observed_min ? min_t_ : 0));

        Diag d_back = std::min<Diag>(get_max_diag(), (observed_max ? max_t_ : static_cast<SOffset>(tsize_)) - (q_e - 1));
        assert(d_back + q_e - 1 <= (observed_max ? max_t_ : static_cast<SOffset>(tsize_)));

        return std::make_pair(d_front, d_back);
    }

    std::vector<std::pair<SOffset, SOffset>> get_overlaps_with_q(SOffset q_b, SOffset q_e) {
        q_b = std::max(q_b, min_q_);
        q_e = std::min(q_e, max_q_ + 1);
        std::vector<std::pair<SOffset, SOffset>> overlaps;
        for (SOffset q = q_b; q < q_e; ++q) {
            if (qpoints_[q] != qpoints_.size()) {
                Diag k = idx_to_diag(qpoints_[q]);
                SOffset t = k + q;
                assert(t >= min_t_);
                assert(t <= max_t_);
                overlaps.emplace_back(q, t);
            }
        }
        return overlaps;
    }

    template <bool observed_min = false, bool observed_max = false>
    std::pair<Diag, Diag> get_diag_range_from_t(SOffset t_b, SOffset t_e) {
        // k = t - q => q = t - k
        // we need q >= q_obs_min and q <= t_obs_max
        // t - k >= q_obs_min => t - q_obs_min >= k
        // t - k <= q_obs_max => t - q_obs_max <= k
        Diag d_front = std::max<Diag>(get_min_diag(), t_e - 1 - (observed_max ? max_q_ : static_cast<SOffset>(qsize_)));
        assert(t_e - d_front - 1 <= (observed_max ? max_q_ : static_cast<SOffset>(qsize_)));

        Diag d_back = std::min<Diag>(get_max_diag(), t_b - (observed_min ? min_q_ : 0));
        assert(t_b - d_back >= (observed_min ? min_q_ : 0));

        return std::make_pair(d_front, d_back);
    }

    std::vector<std::pair<SOffset, SOffset>> get_overlaps_with_t(SOffset t_b, SOffset t_e) {
        t_b = std::max(t_b, min_t_);
        t_e = std::min(t_e, max_t_ + 1);
        std::vector<std::pair<SOffset, SOffset>> overlaps;
        for (SOffset t = t_b; t < t_e; ++t) {
            if (tpoints_[t] != tpoints_.size()) {
                Diag k = idx_to_diag(tpoints_[t]);
                SOffset q = t - k;
                assert(q <= max_q_);
                assert(q >= min_q_);
                overlaps.emplace_back(q, t);
            }
        }
        return overlaps;
    }

    Penalty get_global_min(Diag diag) const { return global_min_[diag_to_idx(diag)]; }

    void extend(const Callback &callback = [](Offset, Offset, Offset) {},
                Penalty min_p = ScoreModel::inf_p) {
        wf_extend(qbegin_, qend_, tbegin_, tend_, p_, table_[p_], mismatch_getter,
            [&](Offset q_mm_b, Offset q_mm_e, Offset t_mm_b, Offset t_mm_e,
                Offset, Offset, Offset, Offset) {
                process(min_p, q_mm_b, q_mm_e, t_mm_b, t_mm_e, callback);
            },
            false
        );
    }

    size_t next() {
        ++p_;
        size_t wf_width = wf_next(qbegin_, qend_, tbegin_, tend_, table_, p_, score_model_,
            [this](Offset q, Offset t, UDiag wf_width) {
                Diag diag = get_diag(q, t);

                if (get_global_min(diag) >= ScoreModel::inf_p)
                    return true;

                const auto &cur_data = diag_data(diag);
                if (cur_data.size() && (cur_data.end() - 1)->q2 > static_cast<SOffset>(q))
                    return true;

                return wf_width >= min_wf_len_
                            && static_cast<Diag>(q + t + max_dist_) < frp_q_ + frp_t_;
            }
        );

        if (!wf_width) {
            ++empty_counter_;
        } else {
            empty_counter_ = 0;
        }

        return wf_width;
    }

    void disable_diag(Diag diag) {
        size_t idx = diag_to_idx(diag);
        global_min_[idx] = ScoreModel::inf_p;
        endpoints_[idx].clear();
    }

  private:
    static constexpr auto mismatch_getter = get_mismatch_it<string_const_iterator,
                                                            string_const_iterator>;

    const ScoreModel &score_model_;
    string_const_iterator qbegin_;
    string_const_iterator qend_;
    const size_t qsize_;

    string_const_iterator tbegin_;
    string_const_iterator tend_;
    const size_t tsize_;

    DPTable table_;
    size_t empty_counter_;
    Penalty p_;

    Diag min_k_;
    Diag max_k_;

    SOffset min_q_;
    SOffset max_q_;
    SOffset min_t_;
    SOffset max_t_;

    SOffset frp_q_;
    SOffset frp_t_;
    Penalty frp_p_;
    Offset frp_ext_;

    std::vector<Container> endpoints_;
    std::vector<size_t> qpoints_;
    std::vector<size_t> tpoints_;
    std::vector<Penalty> global_min_;

    const UDiag min_wf_len_;
    const Offset max_dist_;

    size_t diag_to_idx(Diag diag) const {
        assert(static_cast<Diag>(qsize_) + diag >= 0);
        assert(qsize_ + diag < endpoints_.size());
        return qsize_ + diag;
    }

    Diag idx_to_diag(size_t idx) const {
        assert(idx < endpoints_.size());

        Diag k = static_cast<Diag>(idx) - qsize_;
        assert(idx == diag_to_idx(k));

        return k;
    }

    void process(Penalty min_p = ScoreModel::inf_p,
                 SOffset q = 0,
                 SOffset q_e = 0,
                 SOffset t = 0,
                 SOffset t_e = 0,
                 const Callback &callback = [](Offset, Offset, Offset) {}) {
        assert(q_e <= static_cast<SOffset>(qsize_));
        assert(t_e <= static_cast<SOffset>(tsize_));

        Diag diag = get_diag(q_e, t_e);
        min_k_ = std::min(min_k_, diag);
        max_k_ = std::max(max_k_, diag);

        min_q_ = std::min(min_q_, q);
        max_q_ = std::max(max_q_, q_e);
        min_t_ = std::min(min_t_, t);
        max_t_ = std::max(max_t_, t_e);

        Penalty p = p_;

        // TODO: inv_ext_p != 0 not supported yet
        if constexpr (inv_ext) {
            assert(score_model_.inv_ext_p == 0);
            p += q * score_model_.inv_ext_p;
        }

        if constexpr (inv_open) {
            p += score_model_.inv_open_p;
        }

        // update greatest antidiagonal
        SOffset ext = q_e - q;
        if (q_e + t_e > static_cast<SOffset>(frp_q_ + frp_t_)) {
            frp_q_ = q_e;
            frp_t_ = t_e;
            frp_p_ = p;
            frp_ext_ = q_e - q;
        }

        if (p >= min_p) {
            disable_diag(get_diag(q, t));
            return;
        }

        size_t idx = diag_to_idx(get_diag(q, t));

        // Note: if this diagonal has been disabled, then global_min_[idx] == inf_p, so we skip it
        if (global_min_[idx] >= min_p)
            return;

        auto &bucket = endpoints_[idx];
        assert(bucket.empty() || p >= (bucket.end() - 1)->penalty);
        assert(bucket.empty() || p >= global_min_[idx]);

        if (bucket.empty())
            global_min_[idx] = p;

        bool ret_val = bucket.extend(q, q_e + 1, p);
        assert(!bucket.empty());

        if (!ret_val)
            return;

        for (SOffset i = 0; i <= ext; ++i) {
            if (qpoints_[q + i] == qpoints_.size())
                qpoints_[q + i] = idx;

            if (tpoints_[t + i] == tpoints_.size())
                tpoints_[t + i] = idx;
        }

        callback(q, t, ext);
    }
};

std::tuple<Score, std::string, SOffset, SOffset, SOffset, SOffset,
           Score, std::string, SOffset, SOffset, SOffset, SOffset>
run_alignment(wfa::WFAligner& aligner,
              wfa::WFAligner& aligner_scorer,
              const ScoreModel &score_model,
              std::string_view query_1,
              std::string_view query_rc_1,
              std::string_view target_1,
              std::string_view query_2,
              std::string_view query_rc_2,
              std::string_view target_2,
              SOffset heuristics_length_cutoff = max_offset,
              SOffset min_wavefront_length = max_diag_width,
              SOffset max_distance_threshold = max_offset);