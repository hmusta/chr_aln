#pragma once

#include "helpers.hpp"

#include <algorithm>
#include <functional>
#include <numbers>
#include <numeric>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <omp.h>

#include "cpp_progress_bar/progress_bar.hpp"

using kmer_t = uint64_t;

struct Ranges {
    Ranges(SOffset rb, SOffset re, bool ro, SOffset qb, SOffset qe, bool qo, SOffset = 0, SOffset = 0)
        : rbegin(rb), rend(re), rorientation(ro),
          qbegin(qb), qend(qe), qorientation(qo),
          left_trim(0), right_trim(0) { assert_valid(); }

    static_assert(std::is_same_v<SOffset, std::make_signed_t<SOffset>>);
    SOffset rbegin;
    SOffset rend;
    bool rorientation;
    SOffset qbegin;
    SOffset qend;
    bool qorientation;
    SOffset left_trim;
    SOffset right_trim;

    size_t size() const {
        assert_valid();
        return rend - rbegin;
    }

    Diag diag() const {
        assert_valid();
        return rend - qend;
    }

    void assert_valid(bool allow_empty = true) const {
        assert(!rorientation);
        assert(rbegin >= 0);
        assert(qbegin >= 0);
        assert(rbegin <= rend);
        assert(allow_empty || rbegin < rend);
        assert(qbegin <= qend);
        assert(allow_empty || qbegin < qend);
        assert(rend - rbegin == qend - qbegin);
    }

    void shift_start(size_t target_shift, size_t query_shift) {
        assert_valid();
        rbegin += target_shift;
        rend += target_shift;
        qbegin += query_shift;
        qend += query_shift;
    }

    void trim_prefix(size_t trim, bool allow_empty = true) {
        assert_valid();
        assert(trim <= size());
        if (trim > 0) {
            rbegin += trim;
            qbegin += trim;
            assert_valid(allow_empty);
            left_trim += trim;
        }
    }

    void trim_suffix(size_t trim, bool allow_empty = true) {
        assert_valid();
        assert(trim <= size());
        if (trim > 0) {
            rend -= trim;
            qend -= trim;
            assert_valid(allow_empty);
            right_trim += trim;
        }
    }

    std::string_view get_seq_from_target(std::string_view target) const {
        assert_valid();
        assert(rbegin + size() <= target.size());
        return std::string_view(target.data() + rbegin, size());
    }

    std::string_view get_seq_from_query(std::string_view query,
                                        std::string_view query_rc) const {
        assert_valid();
        assert(qorientation || qbegin + size() <= query.size());
        assert(!qorientation || static_cast<ssize_t>(query_rc.size()) >= qend);
        assert(!qorientation || query_rc.size() - qend + size() <= query_rc.size());
        return std::string_view(!qorientation
                                        ? query.data() + qbegin
                                        : query_rc.data() + query_rc.size() - qend,
                                 size());
    }

    bool check_equal(std::string_view target,
                     std::string_view query,
                     std::string_view query_rc) const {
        assert_valid();
        return get_seq_from_target(target) == get_seq_from_query(query, query_rc);
    }
};

inline bool operator<(const Ranges &a, const Ranges &b) {
    bool in_inv = a.qorientation && b.qorientation;

    return std::tie(a.rend, a.rbegin, !in_inv ? a.qend : b.qbegin, !in_inv ? a.qbegin : b.qend)
         < std::tie(b.rend, b.rbegin, !in_inv ? b.qend : a.qbegin, !in_inv ? b.qbegin : a.qend);
}

inline bool operator==(const Ranges &a, const Ranges &b) {
    return std::tie(a.rbegin, a.rend, a.rorientation, a.qbegin, a.qend, a.qorientation)
            == std::tie(b.rbegin, b.rend, b.rorientation, b.qbegin, b.qend, b.qorientation);
}

inline std::ostream& operator<<(std::ostream &out, const Ranges &a) {
    out << orientation_to_char(a.rorientation) << " "
        << a.rbegin + 1 << "-" << a.rend << "\t"
        << orientation_to_char(a.qorientation) << " "
        << a.qbegin + 1 << "-" << a.qend << "\t"
        << a.size() << "(" << a.left_trim << "," << a.right_trim << ")\t"
        << a.rend - a.qend;

    return out;
}

struct ChainScoreModel {
    ChainScoreModel(size_t min_len,
                    size_t gap_switch,
                    double exp_mismatch_frac_between_mum_bp = 0.01,
                    size_t mx_gap = max_offset)
          : match_s(1),
            c1(std::exp(-exp_mismatch_frac_between_mum_bp * min_len)),
            c2(0.05 * c1),
            max_gap(mx_gap),
            inv_ext_s(0),
            inv_open_s(compute_inv_open_s(gap_switch)) {
        assert(c1 > 0);
        assert(inv_ext_s <= 0);
        assert(inv_open_s <= 0);
    }

    Score get_match_score(SOffset len) const {
        assert(len >= 0);
        return len * match_s;
    }

    Score get_mismatch_score(SOffset len) const {
        assert(len >= 0);

        Score mismatch_score = -std::ceil(c2 * len);
        assert(mismatch_score <= 0);
        assert(len == 0 || mismatch_score < 0);

        return mismatch_score;
    }

    Score get_gap_score(Offset gap_edits) const {
        Score gap_score = gap_edits > 0
            ? -std::ceil(get_gap_cost(gap_edits))
            : 0;
        assert(gap_score <= 0);
        assert(gap_edits == 0 || gap_score < 0);
        return gap_score;
    }

    double get_gap_cost(double gap_edits) const {
        assert(gap_edits >= 1.0);
        return c1 * gap_edits + std::log2(gap_edits);
    }

    Score compute_inv_open_s(size_t gap_switch) const {
        // if gap_switch == max_offset, then our scoring model has linear gap cost,
        // so use the beginning of the curve to estimate the slope
        // otherwise, estimate the slope from the long-gap part of the function
        size_t base_bp = gap_switch < max_offset ? gap_switch : 1;

        // slope = get_gap_cost'(base_bp)
        double inv_slope = c1 + std::numbers::ln2_v<double> / base_bp;

        // get y-intercept
        double inv_open = get_gap_cost(base_bp) - inv_slope * base_bp;

        return -std::ceil(inv_open);
    }

    Score match_s;
    double c1;
    double c2;
    SOffset max_gap;

    Score inv_ext_s;
    Score inv_open_s;
};

inline std::ostream& operator<<(std::ostream &out, const ChainScoreModel &a) {
    out << "chain params: "
        << "match: " << a.match_s << "\t"
        << "c1: " << a.c1 << "\t"
        << "c2 (exp. num mismatches): " << a.c2 << "\t"
        << "max gap: " << a.max_gap << "\t"
        << "inv_open_s: " << a.inv_open_s << "\t"
        << "inv_ext_s: " << a.inv_ext_s;

    return out;
}

void print_chain(std::ostream &out,
                 std::string_view theader,
                 std::string_view target,
                 std::string_view qheader,
                 std::string_view query,
                 const std::vector<Ranges> &chain,
                 Score chain_score);

enum LastOp : uint8_t { CLIPPED, INSERT, DELETE, MATCH, MISMATCH, INVERTED };
static_assert(sizeof(LastOp) == 1);

inline char op_to_char(LastOp op) {
    switch (op) {
        case LastOp::CLIPPED: { return 'S'; } break;
        case LastOp::INSERT: { return 'I'; } break;
        case LastOp::DELETE: { return 'D'; } break;
        case LastOp::MATCH: { return '='; } break;
        case LastOp::MISMATCH: { return 'X'; } break;
        case LastOp::INVERTED: { return 'i'; } break;
    }

    assert(false && "This should not happen");
    return 0;
}

struct ChainTableElem {
    static constexpr size_t npos = std::numeric_limits<size_t>::max() / 2;
    using It = typename std::vector<Ranges>::const_iterator;

    ChainTableElem(Score cs, Score csi, It i,
                   size_t l = ChainTableElem::npos,
                   size_t li = ChainTableElem::npos,
                   size_t si = ChainTableElem::npos,
                   LastOp lo = LastOp::CLIPPED,
                   LastOp loi = LastOp::CLIPPED)
        : chain_score(cs), chain_score_inv(csi), it(i),
          last(l), last_inv(li), start_inv(si),
          last_op(lo), last_op_inv(loi) { assert_valid(); }

    void assert_valid() const {
        assert((last == npos) == (last_op == LastOp::CLIPPED));
        assert((last_inv == npos) == (last_op_inv == LastOp::CLIPPED));
        assert(last == npos || last_op == LastOp::MATCH || last_op == LastOp::INVERTED);
        assert(last_inv == npos || last_op_inv == LastOp::MATCH || last_op_inv == LastOp::INVERTED);
    }

    Score chain_score;
    Score chain_score_inv;
    It it;
    size_t last;
    size_t last_inv;
    size_t start_inv;
    LastOp last_op;
    LastOp last_op_inv;
};

inline bool operator<(const ChainTableElem &a, const ChainTableElem &b) {
    return *a.it < *b.it;
}

std::tuple<bool, std::string_view, std::string_view, SOffset, SOffset, SOffset, SOffset>
extract_gap_seqs_continue(std::string_view target,
                          std::string_view query,
                          std::string_view query_rc,
                          const std::vector<Ranges>& best_chain,
                          size_t i);

template <bool is_close>
inline std::tuple<bool, std::string_view, std::string_view, std::string_view, SOffset, SOffset, SOffset, SOffset>
extract_gap_seqs_switch(std::string_view target,
                        std::string_view query,
                        std::string_view query_rc,
                        const std::vector<Ranges>& ranges,
                        size_t i_prev,
                        size_t i_start,
                        size_t i_end,
                        size_t i_next = ChainTableElem::npos) {
    if (i_next == ChainTableElem::npos)
        i_next = ranges.size();

    assert(i_prev < ranges.size());
    assert(i_start < ranges.size());
    assert(i_end < ranges.size());
    assert(i_next <= ranges.size());
    assert(!is_close || i_next < ranges.size());

    const auto& [rbegin_prev, rend_prev, rorientation_prev, qbegin_prev, qend_prev, qorientation_prev, left_trim_prev,
                 right_trim_prev] = ranges[i_prev];
    const auto& [rbegin_start, rend_start, rorientation_start, qbegin_start, qend_start, qorientation_start, left_trim_start,
                 right_trim_start] = ranges[i_start];
    const auto& [rbegin_end, rend_end, rorientation_end, qbegin_end, qend_end, qorientation_end, left_trim_end,
                 right_trim_end] = ranges[i_end];

    assert(qorientation_prev != qorientation_start);

    const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation, left_trim,
                 right_trim] = ranges[!is_close ? i_start : i_next];

    SOffset mum_length = rend - rbegin;
    SOffset target_begin;
    SOffset target_end;

    if constexpr (!is_close) {
        target_begin = rend_prev;
        target_end = std::max(rbegin, rend_prev);
    } else {
        assert(qorientation_end != qorientation);
        target_begin = std::min(rend_end, rbegin);
        target_end = rbegin;
    }

    assert(target_begin <= target_end);

    std::string_view target_w = target.substr(target_begin, target_end - target_begin);

    std::string_view query_w;
    std::string_view query_rc_w;

    SOffset inv_earliest_qbegin = qend_prev;
    SOffset inv_latest_qbegin = std::max(qbegin_end, qend_prev);
    assert(inv_earliest_qbegin <= inv_latest_qbegin);

    SOffset inv_earliest_qend;
    SOffset inv_latest_qend;
    SOffset query_begin;

    // if there's an overlap, forward takes priority
    if constexpr (!is_close) {
        if (i_next < ranges.size()) {
            const auto& [rbegin_next, rend_next, rorientation_next, qbegin_next, qend_next, qorientation_next, left_trim_next,
                         right_trim_next] = ranges[i_next];
            assert(qorientation_end != qorientation_next);
            // assert(rend_end <= rbegin_next);
            // assert(qend_start <= qbegin_next);
            // assert(rbegin_end < rbegin_next && rend_end < rend_next);
            // assert(qbegin_start < qbegin_next && qend_start < qend_next);
            inv_earliest_qend = std::min(qend_start, qbegin_next);
            assert(inv_latest_qbegin < inv_earliest_qend);

            inv_latest_qend = qbegin_next;
        } else {
            inv_earliest_qend = qend_start;
            assert(inv_latest_qbegin < inv_earliest_qend);

            inv_latest_qend = inv_earliest_qend;
        }
        assert(inv_earliest_qend <= inv_latest_qend);

        query_rc_w = query_rc.substr(query_rc.size() - inv_latest_qend, inv_latest_qend - inv_earliest_qend);

        query_w = query.substr(inv_earliest_qbegin, inv_latest_qbegin - inv_earliest_qbegin);
        query_begin = inv_earliest_qbegin;
    } else {
        const auto& [rbegin_next, rend_next, rorientation_next, qbegin_next, qend_next, qorientation_next, left_trim_next,
                     right_trim_next] = ranges[i_next];
        assert(qorientation_end != qorientation_next);

        inv_earliest_qend = std::min(qend_start, qbegin_next);
        assert(inv_latest_qbegin < inv_earliest_qend);

        inv_latest_qend = qbegin_next;
        assert(inv_earliest_qend <= inv_latest_qend);

        query_rc_w = query.substr(inv_earliest_qend, inv_latest_qend - inv_earliest_qend);

        query_w = query_rc.substr(query_rc.size() - inv_latest_qbegin, inv_latest_qbegin - inv_earliest_qbegin);
        query_begin = inv_earliest_qend;
    }

    SOffset o1 = rend - std::max(rbegin, target_begin);
    SOffset o2 = qend - std::max(qbegin, query_begin);
    SOffset overlap = std::min(o1, o2);

    return std::make_tuple(!qorientation,
                           query_w, query_rc_w,
                           target_w,
                           mum_length,
                           query_begin, target_begin,
                           overlap);
}

std::pair<std::vector<Ranges>, Score> chain_ranges(std::string_view target,
                                                   std::string_view query,
                                                   std::string_view query_rc,
                                                   std::vector<Ranges>& ranges,
                                                   const ScoreModel &score_model,
                                                   const ChainScoreModel &chain_score_model,
                                                   bool chain_inversions = false,
                                                   size_t nthreads = 1);

std::vector<Ranges> rc_ranges(const std::vector<Ranges>& mummer_ranges, SOffset query_size);

void reseed_large_gaps(std::string_view target,
                       std::string_view query,
                       std::string_view query_rc,
                       const ScoreModel &score_model,
                       std::vector<Ranges>& best_chain,
                       std::vector<size_t>& inv_starts,
                       std::vector<size_t>& inv_ends,
                       SOffset k,
                       const std::function<bool(SOffset, SOffset)> &check_if_heuristics,
                       SOffset max_gap = max_offset,
                       double exp_mismatch_frac_between_mum_bp = 0.01,
                       bool check_inversions = false,
                       size_t nthreads = 1);

std::pair<Diag, Diag> compute_diag(const std::string& target,
                                   const std::string& query,
                                   const std::string& query_rc,
                                   const std::vector<Ranges>& best_chain,
                                   const std::vector<size_t>& inv_starts,
                                   const std::vector<size_t>& inv_ends);

std::pair<std::vector<size_t>, std::vector<size_t>>
compute_invs(std::string_view target,
             std::string_view query,
             std::string_view query_rc,
             std::vector<Ranges>& best_chain);