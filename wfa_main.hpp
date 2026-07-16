#pragma once

#include "cigar.hpp"
#include "helpers.hpp"
#include "offset_vector.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <cassert>
#include <cstdint>
#include <cstddef>

#include <bindings/cpp/WFAligner.hpp>

// TODO: replace with unique_ptr for generic object
inline wfa::WFAlignerGapAffine2Pieces
make_aligner(const ScoreModel &score_model,
             wfa::WFAligner::MemoryModel memory_model = wfa::WFAligner::MemoryLow,
             wfa::WFAligner::AlignmentScope scope = wfa::WFAligner::Alignment) {
    wfa::WFAlignerGapAffine2Pieces ret_val(
        score_model.mismatch_p,
        score_model.gap_open_p, score_model.gap_ext_p,
        score_model.gap_open2_p, score_model.gap_ext2_p,
        scope, memory_model
    );

    ret_val.setHeuristicNone();

    return ret_val;
}

inline Penalty get_align_penalty(wfa::WFAligner& aligner_scorer,
                                 std::string_view query,
                                 std::string_view target) {
    aligner_scorer.setHeuristicNone();
    SeqPair view_pair(query, target);
    aligner_scorer.alignEnd2End(match_char, &view_pair, query.size(), target.size());
    assert(aligner_scorer.getAlignmentStatus() == wfa::WFAligner::StatusAlgCompleted);
    assert(aligner_scorer.getAlignmentScore() <= 0);
    return -aligner_scorer.getAlignmentScore();
}

using Element = std::tuple<Offset, // M
                           Offset, // I
                           Offset, // D
                           Offset, // max extension
                           Offset, // insert length
                           Offset // deletion length
                           >;


constexpr Element def_elem(npos, npos, npos, 0, 0, 0);

using Elements = std::pair<Diag, Element>;
using Wavefront = std::vector<Elements>;

using DPTable = OffsetVector<Wavefront>;

inline void check_lengths(const wfa::WFAligner& aligner,
                          Penalty max_pen,
                          uint64_t total_length) {
    using wfa_score_t = std::decay_t<decltype(aligner.getAlignmentScore())>;

    if (static_cast<uint64_t>(max_pen) * total_length
            > static_cast<uint64_t>(std::numeric_limits<wfa_score_t>::max())) {
        std::cerr << "length(Q) + length(T) must be below "
                  << static_cast<uint64_t>(std::numeric_limits<wfa_score_t>::max()) / max_pen
                  << std::endl;
        throw std::runtime_error("Sequences are too long");
    }
}

inline void check_entry(Offset qlen, Offset tlen, const Elements &entry) {
    #ifndef NDEBUG
    const auto& diag = entry.first;
    const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = entry.second;
    assert(q_mm_i != npos);
    assert(q_ext <= q_mm_i);
    assert(q_mm_i <= qlen);

    static_assert(static_cast<SOffset>(npos) == -1);
    assert(q_mm_i - q_ext
        == static_cast<Offset>(std::max({ static_cast<SOffset>(q_mm_i - q_ext),
                                            static_cast<SOffset>(q_ins_i),
                                            static_cast<SOffset>(q_del_i) })));

    assert(static_cast<SOffset>(q_mm_i) + diag >= 0);
    assert(q_mm_i + diag <= tlen);
    assert(q_ext <= q_mm_i + diag);

    if (q_ins_i != npos) {
        assert(nins);
        assert(q_ins_i <= qlen);
        assert(static_cast<SOffset>(q_ins_i) + diag >= 0);
        assert(q_ins_i + diag <= tlen);
    }

    if (q_del_i != npos) {
        assert(ndel);
        assert(q_del_i <= qlen);
        assert(static_cast<SOffset>(q_del_i) + diag >= 0);
        assert(q_del_i + diag <= tlen);
    }
    #endif
}

using WFElemCallback = std::function<void(Offset, Offset, Offset, Offset, Offset, Offset, Offset, Offset)>;

template <class MismatchGetter, class T, class U>
inline bool wf_extend(
        T q_begin,
        T q_end,
        U t_begin,
        U t_end,
        Penalty p,
        Wavefront& wf,
        const MismatchGetter& get_mismatch_it,
        const WFElemCallback& callback = [](Offset, Offset, Offset, Offset, Offset, Offset) {},
        bool halt_if_found = true,
        bool callback_after_extend = false) {
    bool found = false;
    Offset qlen = std::distance(q_begin, q_end);
    Offset tlen = std::distance(t_begin, t_end);
    for (auto& entry : wf) {
        check_entry(qlen, tlen, entry);
        auto& [diag, offsets] = entry;
        auto& [q_mm, q_ins, q_del, q_ext, nins, ndel] = offsets;
        assert(q_ext == 0);

        Offset t_mm = q_mm + diag;

        Offset old_q_mm = q_mm;
        Offset old_t_mm = t_mm;
        q_mm = get_mismatch_it(q_begin + old_q_mm, q_end, t_begin + old_t_mm, t_end) - q_begin;
        assert(q_mm >= old_q_mm);
        q_ext = q_mm - old_q_mm;
        t_mm += q_ext;

        check_entry(qlen, tlen, entry);

        if (!callback_after_extend)
            callback(old_q_mm, q_mm, old_t_mm, t_mm, q_ins, q_del, nins, ndel);

        found |= (q_mm == qlen && t_mm == tlen);

        if (found && halt_if_found) {
            if (callback_after_extend)
                callback(old_q_mm, q_mm, old_t_mm, t_mm, q_ins, q_del, nins, ndel);

            return true;
        }
    }

    if (callback_after_extend) {
        for (const auto& [diag, offsets] : wf) {
            const auto& [q_mm, q_ins, q_del, q_ext, nins, ndel] = offsets;
            Offset t_mm = q_mm + diag;
            assert(t_mm <= tlen);
            callback(q_mm - q_ext, q_mm, t_mm - q_ext, t_mm, q_ins, q_del, nins, ndel);
        }
    }

    return found;
}

static const Wavefront dummy_wavefront;

using DiagSkipper = std::function<bool(Diag, UDiag, UDiag)>;

template <class T, class U>
inline size_t wf_next(T q_begin,
                      T q_end,
                      U t_begin,
                      U t_end,
                      DPTable& table,
                      Penalty p,
                      const ScoreModel &score_model,
                      const DiagSkipper &skip_diag = [](Offset, Offset, UDiag) { return false; }) {
    Offset qlen = std::distance(q_begin, q_end);
    Offset tlen = std::distance(t_begin, t_end);
    assert(p <= table.size());

    if (p == table.size())
        table.emplace_back();

    assert(p >= table.offset());
    assert(table.size() == table.alloc_size() + table.offset());

    if (p - table.offset() > score_model.max_pen) {
        table.dealloc_front(std::min(p - table.offset() - score_model.max_pen, table.alloc_size() - 1));
        assert(table.alloc_size() == score_model.max_pen + 1);
    }

    assert(table.alloc_size());
    assert(p < table.size());
    assert(p >= table.offset());

    auto mt = dummy_wavefront.cbegin();
    auto mt_end = dummy_wavefront.cend();

    auto eit = dummy_wavefront.cbegin();
    auto eit_end = dummy_wavefront.cend();

    auto edt = dummy_wavefront.cbegin();
    auto edt_end = dummy_wavefront.cend();

    auto oit = dummy_wavefront.cbegin();
    auto oit_end = dummy_wavefront.cend();

    auto odt = dummy_wavefront.cbegin();
    auto odt_end = dummy_wavefront.cend();

    auto e2it = dummy_wavefront.cbegin();
    auto e2it_end = dummy_wavefront.cend();

    auto e2dt = dummy_wavefront.cbegin();
    auto e2dt_end = dummy_wavefront.cend();

    auto o2it = dummy_wavefront.cbegin();
    auto o2it_end = dummy_wavefront.cend();

    auto o2dt = dummy_wavefront.cbegin();
    auto o2dt_end = dummy_wavefront.cend();

    auto& wf = table[p];

    auto it = wf.begin();
    auto it_end = wf.end();

    Diag obs_diag_min = max_diag;
    Diag obs_diag_max = min_diag;

    Diag mn = -static_cast<Diag>(qlen);
    Diag mx = tlen;

    auto check_diags = [mn,mx](Diag d_front, Diag d_back) -> bool {
        return d_front <= d_back && d_front >= mn && d_back <= mx;
    };

    auto update_diag = [&obs_diag_min,&obs_diag_max,&check_diags,mn,mx]<Diag delta>(typename Wavefront::const_iterator &d_begin,
                                                                                    typename Wavefront::const_iterator &d_end) {
        assert(d_begin != d_end);
        assert(check_diags(d_begin->first, (d_end - 1)->first));
        static_assert(delta >= -1);
        static_assert(delta <= 1);

        if constexpr (delta == -1) {
            if (d_begin->first - 1 < mn)
                ++d_begin;

            if (d_begin != d_end) {
                assert((d_end - 1)->first - 1 >= mn);
                obs_diag_min = std::min(obs_diag_min, d_begin->first - 1);
                obs_diag_max = std::max(obs_diag_max, (d_end - 1)->first - 1);
                assert(check_diags(obs_diag_min, obs_diag_max));
            }
        } else if constexpr (delta == 1) {
            if (d_begin != d_end && (d_end - 1)->first + 1 > mx)
                --d_end;

            if (d_begin != d_end) {
                assert(d_begin->first + 1 <= mx);
                obs_diag_min = std::min(obs_diag_min, d_begin->first + 1);
                obs_diag_max = std::max(obs_diag_max, (d_end - 1)->first + 1);
                assert(check_diags(obs_diag_min, obs_diag_max));
            }
        } else {
            static_assert(delta == 0 && "Only values of -1, 0, and 1 are allowed");
            obs_diag_min = std::min(obs_diag_min, d_begin->first);
            obs_diag_max = std::max(obs_diag_max, (d_end - 1)->first);
        }
    };

    if (wf.size()) {
        assert(p == 0);
        assert(wf.size() == 1);
        obs_diag_min = it->first;
        obs_diag_max = (it_end - 1)->first;
        assert(check_diags(obs_diag_min, obs_diag_max));
    }

    Penalty last_mismatch_p = p - score_model.mismatch_p;
    assert(p < score_model.mismatch_p || last_mismatch_p < table.size());
    if (p >= score_model.mismatch_p && table[last_mismatch_p].size()) {
        const auto &last_mismatch_table = table[last_mismatch_p];
        mt = last_mismatch_table.begin();
        mt_end = last_mismatch_table.end();
        update_diag.template operator()<0>(mt, mt_end);
    }

    Penalty last_ext_p = p - score_model.gap_ext_p;
    assert(p < score_model.gap_ext_p || last_ext_p < table.size());
    if (p >= score_model.gap_ext_p && table[last_ext_p].size()) {
        const auto &last_ext_table = table[last_ext_p];
        eit = last_ext_table.begin();
        eit_end = last_ext_table.end();
        update_diag.template operator()<-1>(eit, eit_end);

        edt = last_ext_table.begin();
        edt_end = last_ext_table.end();
        update_diag.template operator()<1>(edt, edt_end);
    }

    Penalty last_open_p = p - score_model.get_gap_penalty(1);
    assert(p < score_model.gap_open_p + score_model.gap_ext_p || last_open_p < table.size());
    if (p >= score_model.gap_open_p + score_model.gap_ext_p && table[last_open_p].size()) {
        const auto &last_open_table = table[last_open_p];
        oit = last_open_table.begin();
        oit_end = last_open_table.end();
        update_diag.template operator()<-1>(oit, oit_end);

        odt = last_open_table.begin();
        odt_end = last_open_table.end();
        update_diag.template operator()<1>(odt, odt_end);
    }

    if (score_model.penalty_diff < ScoreModel::inf_p) {
        Penalty last_ext2_p = p - score_model.gap_ext2_p;
        assert(p < score_model.gap_ext2_p || last_ext2_p < table.size());
        if (p >= score_model.gap_ext2_p && table[last_ext2_p].size()) {
            const auto &last_ext2_table = table[last_ext2_p];
            e2it = last_ext2_table.begin();
            e2it_end = last_ext2_table.end();
            update_diag.template operator()<-1>(e2it, e2it_end);

            e2dt = last_ext2_table.begin();
            e2dt_end = last_ext2_table.end();
            update_diag.template operator()<1>(e2dt, e2dt_end);
        }

        Penalty last_gap_diff_p = p - score_model.penalty_diff;
        assert(p < score_model.penalty_diff || last_gap_diff_p < table.size());
        if (p >= score_model.penalty_diff && table[last_gap_diff_p].size()) {
            const auto &last_gap_diff_table = table[last_gap_diff_p];
            o2it = last_gap_diff_table.begin();
            o2it_end = last_gap_diff_table.end();
            update_diag.template operator()<-1>(o2it, o2it_end);

            o2dt = last_gap_diff_table.begin();
            o2dt_end = last_gap_diff_table.end();
            update_diag.template operator()<1>(o2dt, o2dt_end);
        }
    }

    if (obs_diag_min > obs_diag_max)
        return 0;

    assert(check_diags(obs_diag_min, obs_diag_max));

    UDiag wf_width = obs_diag_max - obs_diag_min + 1;
    Wavefront merged_wavefront;
    merged_wavefront.reserve(obs_diag_max - obs_diag_min + 1);



    for (Diag diag = obs_diag_min; diag <= obs_diag_max; ++diag) {
        assert(merged_wavefront.empty() || merged_wavefront.back().first < diag);
        if (it != it_end && it->first == diag) {
            assert(p == 0);
            check_entry(qlen, tlen, *it);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = it->second;
            if (!skip_diag(q_mm_i, diag + q_mm_i, wf_width))
                merged_wavefront.emplace_back(*it);

            ++it;
            assert(it == it_end);
        }

        // mismatch
        if (mt != mt_end && mt->first == diag) {
            check_entry(qlen, tlen, *mt);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = mt->second;
            auto t_mm_i = q_mm_i + diag;

            if (q_mm_i < qlen && t_mm_i < tlen && !skip_diag(q_mm_i + 1, t_mm_i + 1, wf_width)) {
                assert(*(q_begin + q_mm_i) != *(t_begin + t_mm_i));
                Offset reach_q = q_mm_i + 1;
                assert(reach_q != npos);
                assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);
                if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                    std::get<0>(merged_wavefront.emplace_back(diag, def_elem).second)
                            = reach_q;
                } else {
                    auto& val = std::get<0>(merged_wavefront.back().second);
                    if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                        val = reach_q;
                    }
                }
                assert(std::get<0>(merged_wavefront.back().second) != npos);
            }
            ++mt;
        }

        // insert extend
        if (eit != eit_end && eit->first - 1 == diag) {
            check_entry(qlen, tlen, *eit);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = eit->second;

            if (q_ins_i != npos && static_cast<SOffset>(nins) + 1 < score_model.gap_switch) {
                Offset reach_q = q_ins_i + 1;

                auto t_ins_i = reach_q + diag;
                assert(t_ins_i == q_ins_i + eit->first);

                if (reach_q <= qlen && !skip_diag(reach_q, t_ins_i, wf_width)) {
                    assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                    if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                        std::get<1>(merged_wavefront.emplace_back(diag, def_elem).second)
                                = reach_q;
                        std::get<4>(merged_wavefront.back().second) = nins + 1;
                    } else {
                        auto& val = std::get<1>(merged_wavefront.back().second);
                        if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                            val = reach_q;
                            std::get<4>(merged_wavefront.back().second) = nins + 1;
                        }
                    }
                    assert(std::get<1>(merged_wavefront.back().second) != npos);
                }
            }

            ++eit;
        }

        // insert extend 2
        if (e2it != e2it_end && e2it->first - 1 == diag) {
            check_entry(qlen, tlen, *e2it);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = e2it->second;

            if (q_ins_i != npos && static_cast<SOffset>(nins) >= score_model.gap_switch) {
                Offset reach_q = q_ins_i + 1;

                Offset t_ins_i = reach_q + diag;
                assert(t_ins_i == q_ins_i + e2it->first);

                if (reach_q <= qlen && !skip_diag(reach_q, t_ins_i, wf_width)) {
                    assert(reach_q != npos);
                    assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                    if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                        std::get<1>(merged_wavefront.emplace_back(diag, def_elem).second)
                                = reach_q;
                        std::get<4>(merged_wavefront.back().second) = nins + 1;
                    } else {
                        auto& val = std::get<1>(merged_wavefront.back().second);
                        if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                            val = reach_q;
                            std::get<4>(merged_wavefront.back().second) = nins + 1;
                        }
                    }
                    assert(std::get<1>(merged_wavefront.back().second) != npos);
                }
            }

            ++e2it;
        }

        // delete extend
        if (edt != edt_end && edt->first + 1 == diag) {
            check_entry(qlen, tlen, *edt);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = edt->second;

            if (q_del_i != npos && static_cast<SOffset>(ndel) + 1 < score_model.gap_switch) {
                auto t_del_i = q_del_i + diag - 1;
                assert(t_del_i == q_del_i + edt->first);

                if (t_del_i + 1 <= tlen && !skip_diag(q_del_i, t_del_i + 1, wf_width)) {
                    assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                    Offset reach_q = q_del_i;
                    if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                        std::get<2>(merged_wavefront.emplace_back(diag, def_elem).second)
                                = reach_q;
                        std::get<5>(merged_wavefront.back().second) = ndel + 1;
                    } else {
                        auto& val = std::get<2>(merged_wavefront.back().second);
                        if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                            val = reach_q;
                            std::get<5>(merged_wavefront.back().second) = ndel + 1;
                        }
                    }
                    assert(std::get<2>(merged_wavefront.back().second) != npos);
                }
            }

            ++edt;
        }

        // delete extend 2
        if (e2dt != e2dt_end && e2dt->first + 1 == diag) {
            check_entry(qlen, tlen, *e2dt);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = e2dt->second;

            if (q_del_i != npos && static_cast<SOffset>(ndel) >= score_model.gap_switch) {
                auto t_del_i = q_del_i + diag - 1;
                assert(t_del_i == q_del_i + e2dt->first);

                if (t_del_i + 1 <= tlen && !skip_diag(q_del_i, t_del_i + 1, wf_width)) {
                    assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                    Offset reach_q = q_del_i;
                    if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                        std::get<2>(merged_wavefront.emplace_back(diag, def_elem).second)
                                = reach_q;
                        std::get<5>(merged_wavefront.back().second) = ndel + 1;
                    } else {
                        auto& val = std::get<2>(merged_wavefront.back().second);
                        if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                            val = reach_q;
                            std::get<5>(merged_wavefront.back().second) = ndel + 1;
                        }
                    }
                    assert(std::get<2>(merged_wavefront.back().second) != npos);
                }
            }

            ++e2dt;
        }

        // insert open
        if (oit != oit_end && oit->first - 1 == diag) {
            check_entry(qlen, tlen, *oit);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = oit->second;

            Offset reach_q = q_mm_i + 1;

            Offset t_mm_i = reach_q + diag;
            assert(t_mm_i == q_mm_i + oit->first);

            if (reach_q <= qlen && !skip_diag(reach_q, t_mm_i, wf_width)) {
                assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                    std::get<1>(merged_wavefront.emplace_back(diag, def_elem).second)
                            = reach_q;
                    std::get<4>(merged_wavefront.back().second) = 1;
                } else {
                    auto& val = std::get<1>(merged_wavefront.back().second);
                    if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                        val = reach_q;
                        std::get<4>(merged_wavefront.back().second) = 1;
                    }
                }
                assert(std::get<1>(merged_wavefront.back().second) != npos);
            }

            ++oit;
        }

        // insert open 2
        if (o2it != o2it_end && o2it->first - 1 == diag) {
            check_entry(qlen, tlen, *o2it);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = o2it->second;

            if (q_ins_i != npos && static_cast<SOffset>(nins) + 1 == score_model.gap_switch) {
                Offset reach_q = q_ins_i + 1;

                Offset t_ins_i = reach_q + diag;
                assert(t_ins_i == q_ins_i + o2it->first);

                if (reach_q <= qlen && !skip_diag(reach_q, t_ins_i, wf_width)) {
                    assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                    if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                        std::get<1>(merged_wavefront.emplace_back(diag, def_elem).second)
                                = reach_q;
                        std::get<4>(merged_wavefront.back().second) = nins + 1;
                    } else {
                        auto& val = std::get<1>(merged_wavefront.back().second);
                        if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                            val = reach_q;
                            std::get<4>(merged_wavefront.back().second) = nins + 1;
                        }
                    }
                    assert(std::get<1>(merged_wavefront.back().second) != npos);
                }
            }

            ++o2it;
        }

        // delete open
        if (odt != odt_end && odt->first + 1 == diag) {
            check_entry(qlen, tlen, *odt);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = odt->second;

            auto t_mm_i = q_mm_i + diag - 1;
            assert(t_mm_i == q_mm_i + odt->first);

            if (t_mm_i + 1 <= tlen && !skip_diag(q_mm_i, t_mm_i + 1, wf_width)) {
                assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                Offset reach_q = q_mm_i;
                if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                    std::get<2>(merged_wavefront.emplace_back(diag, def_elem).second)
                            = reach_q;
                    std::get<5>(merged_wavefront.back().second) = 1;
                } else {
                    auto& val = std::get<2>(merged_wavefront.back().second);
                    if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                        val = reach_q;
                        std::get<5>(merged_wavefront.back().second) = 1;
                    }
                }
                assert(std::get<2>(merged_wavefront.back().second) != npos);
            }

            ++odt;
        }

        // delete open 2
        if (o2dt != o2dt_end && o2dt->first + 1 == diag) {
            check_entry(qlen, tlen, *o2dt);
            const auto& [q_mm_i, q_ins_i, q_del_i, q_ext, nins, ndel] = o2dt->second;

            if (q_del_i != npos && static_cast<SOffset>(ndel) + 1 == score_model.gap_switch) {
                auto t_del_i = q_del_i + diag - 1;
                assert(t_del_i == q_del_i + o2dt->first);

                if (t_del_i + 1 <= tlen && !skip_diag(q_del_i, t_del_i + 1, wf_width)) {
                    assert(merged_wavefront.empty() || merged_wavefront.back().first <= diag);

                    Offset reach_q = q_del_i;
                    if (merged_wavefront.empty() || merged_wavefront.back().first < diag) {
                        std::get<2>(merged_wavefront.emplace_back(diag, def_elem).second)
                                = reach_q;
                        std::get<5>(merged_wavefront.back().second) = ndel + 1;
                    } else {
                        auto& val = std::get<2>(merged_wavefront.back().second);
                        if (static_cast<SOffset>(reach_q) > static_cast<SOffset>(val)) {
                            val = reach_q;
                            std::get<5>(merged_wavefront.back().second) = ndel + 1;
                        }
                    }
                    assert(std::get<2>(merged_wavefront.back().second) != npos);
                }
            }

            ++o2dt;
        }
    }

    assert(it == it_end);
    assert(mt == mt_end);
    assert(eit == eit_end);
    assert(edt == edt_end);
    assert(e2it == e2it_end);
    assert(e2dt == e2dt_end);
    assert(oit == oit_end);
    assert(odt == odt_end);
    assert(o2it == o2it_end);
    assert(o2dt == o2dt_end);

    std::swap(merged_wavefront, wf);
    assert(std::is_sorted(wf.begin(), wf.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; }));

    for (auto& entry : wf) {
        auto& [diag, offsets] = entry;
        auto& [q_mm, q_ins, q_del, q_ext, nins, ndel] = offsets;
        assert(q_ext == 0);

        static_assert(static_cast<SOffset>(npos) == -1);
        q_mm = static_cast<Offset>(
                std::max({ static_cast<SOffset>(q_mm), static_cast<SOffset>(q_ins),
                           static_cast<SOffset>(q_del) }));
        check_entry(qlen, tlen, entry);
    }

    return wf.size();
}

inline std::vector<std::pair<SOffset, SOffset>> get_ns(std::string_view seq) {
    std::vector<std::pair<SOffset, SOffset>> n_coords;
    for (auto it = seq.begin(); it != seq.end(); ++it) {
        if (*it == 'N') {
            size_t begin = it - seq.begin();
            auto end = seq.find_first_not_of('N', begin);
            if (end == std::string_view::npos) {
                end = seq.size();
            }
            assert(end > begin);

            assert(end == seq.size() || seq[end] != 'N');
            if (end - begin >= 10000) {
                std::cout << "Found gap of length " << end - begin
                          << " in sequence of length " << seq.size() << "\n";
                n_coords.emplace_back(begin, end);
            }

            it = seq.begin() + end - 1;
        }
    }

    return n_coords;
}

inline std::pair<std::string, Score>
run_alignment_gaps(wfa::WFAligner& aligner,
                                      std::string_view query,
                                      std::string_view target,
                                      const ScoreModel &score_model) {
    check_lengths(aligner, score_model.max_pen, query.size() + target.size());
    std::vector<std::pair<SOffset, SOffset>> query_n_coords = get_ns(query);
    std::vector<std::pair<SOffset, SOffset>> target_n_coords = get_ns(target);

    if (query_n_coords.size() != 1 && target_n_coords.size() != 1) {
        std::cout << "Defaulting to WFA2-lib\n";
        SeqPair view_pair(query, target);
        aligner.alignEnd2End(match_char, &view_pair, query.size(), target.size());
        std::string cigar = cigar_fix_n(aligner.getCIGAR(true), target, query);
        Score aln_score = score_cigar(cigar, view_pair, score_model);
        return std::make_pair(cigar, aln_score);
    }

    std::string cigar;

    SOffset query_gap_begin = query_n_coords.size() ? query_n_coords[0].first : query.size();
    SOffset query_gap_end = query_n_coords.size() ? query_n_coords[0].second : query.size();
    SOffset query_gap_length = query_gap_end - query_gap_begin;

    SOffset target_gap_begin
            = target_n_coords.size() ? target_n_coords[0].first : target.size();
    SOffset target_gap_end
            = target_n_coords.size() ? target_n_coords[0].second : target.size();
    SOffset target_gap_length = target_gap_end - target_gap_begin;

    SOffset first_aln_begin = std::min(query_gap_begin, target_gap_begin);
    {
        SeqPair view_pair(std::string_view(query.data(), first_aln_begin),
                          std::string_view(target.data(), first_aln_begin));
        std::cout << "Aligning first " << first_aln_begin << "\n";
        aligner.alignEnd2End(match_char, &view_pair, view_pair.first.size(),
                             view_pair.second.size());
        cigar = cigar_fix_n(aligner.getCIGAR(true), view_pair.second, view_pair.first);
        #ifndef NDEBUG
        Score score = score_cigar(cigar, view_pair, score_model);
        std::cout << "\twith score " << score << "\n";
        #endif
    }

    SOffset qi = first_aln_begin;
    SOffset ti = first_aln_begin;

    SOffset taken_gaps = 0;
    SOffset taken_gaps_r = 0;
    SOffset taken_gaps_q = 0;

    while (taken_gaps_r < target_gap_length || taken_gaps_q < query_gap_length) {
        bool found = false;
        if (qi >= query_gap_begin && qi < query_gap_end) {
            SOffset query_gap_left = query_gap_end - qi;
            assert(query_gap_left > 0);
            if (ti < static_cast<SOffset>(target.size())) {
                SOffset target_left = std::min<SOffset>(query_gap_left, target.size() - ti);
                taken_gaps += target_left;
                qi += target_left;
                taken_gaps_q += target_left;

                if (ti + target_left >= target_gap_end) {
                    taken_gaps_r = target_gap_length;
                } else if (ti <= target_gap_begin && ti + target_left > target_gap_begin) {
                    taken_gaps_r += ti + target_left - target_gap_begin;
                } else if (ti > target_gap_begin && ti + target_left <= target_gap_end) {
                    taken_gaps_r += target_left;
                }

                assert(taken_gaps_r <= target_gap_length);
                ti += target_left;
                found = true;
            }
        } else if (ti >= target_gap_begin && ti < target_gap_end) {
            SOffset target_gap_left = target_gap_end - ti;
            assert(target_gap_length > 0);
            if (qi < static_cast<SOffset>(query.size())) {
                SOffset query_left = std::min<SOffset>(target_gap_left, query.size() - qi);
                taken_gaps += query_left;
                ti += query_left;
                taken_gaps_r += query_left;

                if (qi + query_left >= query_gap_end) {
                    taken_gaps_q = query_gap_length;
                } else if (qi <= query_gap_begin && ti + query_left > query_gap_begin) {
                    taken_gaps_q += ti + query_left - query_gap_begin;
                } else if (qi > query_gap_begin && ti + query_left <= query_gap_end) {
                    taken_gaps_q += query_left;
                }

                assert(taken_gaps_q <= query_gap_length);
                qi += query_left;
                found = true;
            }
        }

        if (!found)
            break;
    }

    if (taken_gaps > 0) {
        assert(ti == first_aln_begin + taken_gaps);
        assert(qi == first_aln_begin + taken_gaps);
        #ifndef NDEBUG
        for (int64_t i = 0; i < taken_gaps; ++i) {
            assert(query[first_aln_begin + i] == 'N' || target[first_aln_begin + 1] == 'N');
        }
        #endif
        std::cout << "Then gap of length " << taken_gaps << "\n";
        cigar += std::to_string(taken_gaps) + "M";

        #ifndef NDEBUG
        Score score = score_cigar(cigar, SeqPair(query.substr(0, qi), target.substr(0, ti)), score_model);
        std::cout << "Score so far score " << score << "\n";
        #endif
    }

    std::cout << "Took gaps " << taken_gaps_r << " / " << target_gap_length << "\t"
              << taken_gaps_q << " / " << query_gap_length << "\n";

    SOffset query_left = query.size() - qi;
    SOffset target_left = target.size() - ti;

    if (query_left == 0 && target_left > 0) {
        std::cout << "Final insertion of length " << target_left << "\n";
        cigar += std::to_string(target_left) + "I";
        #ifndef NDEBUG
        Score score = score_cigar(cigar, SeqPair(query, target), score_model);
        std::cout << "Score so far score " << score << "\n";
        #endif
    } else if (query_left > 0 && target_left == 0) {
        std::cout << "Final deletion of length " << query_left << "\n";
        cigar += std::to_string(query_left) + "D";
        #ifndef NDEBUG
        Score score = score_cigar(cigar, SeqPair(query, target), score_model);
        std::cout << "Score so far score " << score << "\n";
        #endif
    } else if (query_left > 0 && target_left > 0) {
        std::cout << "Final alignment of lengths " << target_left << " " << query_left
                  << "\n";
        SeqPair view_pair(std::string_view(query.data() + qi, query_left),
                          std::string_view(target.data() + ti, target_left));
        aligner.alignEnd2End(match_char, &view_pair, view_pair.first.size(),
                             view_pair.second.size());
        std::string last_cigar = cigar_fix_n(aligner.getCIGAR(true), view_pair.second, view_pair.first);
        #ifndef NDEBUG
        Score score = score_cigar(last_cigar, view_pair, score_model);
        std::cout << "\twith score " << score << "\n";
        #endif
        cigar += std::move(last_cigar);
    }

    Score aln_score = score_cigar(cigar, SeqPair(query, target), score_model);

    return std::make_pair(std::move(cigar), aln_score);
}