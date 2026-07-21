#include "wfa_switch.hpp"

#include "repeat_aligner.hpp"

struct Breakpoints {
    Penalty min_p;
    Offset t1;
    Offset t2;
    Offset q1;
    Offset q2;
    Offset t1_left_gap = 0;
    Offset t1_right_gap = 0;
    Offset t2_left_gap = 0;
    Offset t2_right_gap = 0;
    Offset q1_left_gap = 0;
    Offset q1_right_gap = 0;
    Offset q2_left_gap = 0;
    Offset q2_right_gap = 0;

    std::string_view get_t1_fw(std::string_view target_1) const {
        return target_1.substr(0, t1);
    }
    std::string_view get_t1_bw(std::string_view target_1) const {
        return target_1.substr(t1 + t1_left_gap + t1_right_gap);
    }

    std::string_view get_q1_fw(std::string_view query_1) const {
        return query_1.substr(0, q1);
    }

    std::string_view get_q1_bw(std::string_view query_rc_1) const {
        return query_rc_1.substr(q2 + q2_left_gap + q2_right_gap);
    }

    std::string_view get_t2_fw(std::string_view target_2) const {
        return target_2.substr(0, t2);
    }
    std::string_view get_t2_bw(std::string_view target_2) const {
        return target_2.substr(t2 + t2_left_gap + t2_right_gap);
    }

    std::string_view get_q2_fw(std::string_view query_2) const {
        return query_2.substr(0, query_2.size() - q1 - q1_left_gap - q1_right_gap);
    }

    std::string_view get_q2_bw(std::string_view query_rc_2) const {
        return query_rc_2.substr(query_rc_2.size() - q2);
    }

    static std::tuple<Score, std::string, Offset, Offset>
    align_segment(wfa::WFAligner& aligner,
                  const ScoreModel &score_model,
                  std::string_view query,
                  std::string_view target,
                  bool penalty_to_score = true,
                  SOffset heuristics_length_cutoff = max_offset,
                  Diag min_k = min_diag,
                  Diag max_k = max_diag) {
        SeqPair view_pair = std::make_pair(query, target);
        Score n_penalty = INT32_MIN;
        std::string cigar;
        if (query.size() && target.size()) {
            if (static_cast<SOffset>(query.size() + target.size()) >= heuristics_length_cutoff) {
                cigar = repeat_aligner(std::string(query), std::string(target));
            } else {
                if (min_k == min_diag && max_k == max_diag) {
                    aligner.setHeuristicNone();
                } else {
                    aligner.setHeuristicBandedStatic(-max_k - 1, -min_k + 1);
                }
                aligner.alignEnd2End(match_char, &view_pair, query.size(), target.size());
                assert(aligner.getAlignmentStatus() == wfa::WFAligner::StatusAlgCompleted);
                n_penalty = aligner.getAlignmentScore();
                cigar = cigar_fix_n(aligner.getCIGAR(true), target, query);
                // if (get_cigar) {
                // } else {
                //     cigar = std::to_string(target.size()) + "I" + std::to_string(query.size()) + "D";
                // }
            }
        } else if (query.size()) {
            cigar = std::to_string(query.size()) + "D";
        } else if (target.size()) {
            cigar = std::to_string(target.size()) + "I";
        } else {
            n_penalty = 0;
        }

        assert((query.empty() && target.empty()) || cigar.size());
        assert(check_cigar_seq_lengths(cigar, target.size(), query.size()));
        Score score;
        if (n_penalty != INT32_MIN) {
            score = penalty_to_score
                ? score_model.penalty_to_score(-n_penalty, query.size(), target.size())
                : -n_penalty;
        } else {
            score = score_cigar(cigar, view_pair, score_model, !penalty_to_score);
        }

        return std::make_tuple(score, std::move(cigar), target.size(), query.size());
    }

    bool all_on_right(Offset t1_max, Offset q1_max, Offset t2_max, Offset q2_max) const {
        return q1 + q1_left_gap == 0 && t1 + t1_left_gap == t1_max
                        && q2 == q2_max && (t2_max > 0 || q1_max > 0 || q2_max > 0);
    }

    std::tuple<Score, std::string, Offset, Offset>
    align_all_right(wfa::WFAligner& aligner,
                    const ScoreModel &score_model,
                    std::string_view target_2,
                    std::string_view query_2,
                    std::string_view query_rc_2,
                    bool get_cigar,
                    SOffset heuristics_length_cutoff,
                    bool penalty_to_score) const {
        if (!t2_left_gap && !t2_right_gap && !q1_right_gap && !q2_right_gap) {
            std::string query_2_cat(query_2);
            query_2_cat += query_rc_2;
            auto [p, cigar, r_consumed, q_consumed] = align_segment(
                aligner, score_model, query_2_cat, target_2,
                penalty_to_score,
                heuristics_length_cutoff
            );

            return std::make_tuple(
                p + (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * query_2.size(),
                std::move(cigar),
                r_consumed,
                q_consumed
            );
        } else {
            auto [p2f, cigar2f, rcons2f, qcons2f] = align_segment(
                aligner, score_model, get_q2_fw(query_2), get_t2_fw(target_2),
                penalty_to_score,
                heuristics_length_cutoff
            );
            if (t2_left_gap + t2_right_gap) {
                cigar2f += std::to_string(t2_left_gap + t2_right_gap) + "I";
                p2f += penalty_to_score
                    ? score_model.get_gap_score(t2_left_gap + t2_right_gap)
                    : score_model.get_gap_penalty(t2_left_gap + t2_right_gap);
                rcons2f += t2_left_gap + t2_right_gap;
            }
            if (q1_right_gap + q2_right_gap) {
                cigar2f += std::to_string(q1_right_gap + q2_right_gap) + "D";
                if (penalty_to_score) {
                    p2f += score_model.get_gap_score(q1_right_gap + q2_right_gap)
                        + score_model.inv_ext_s * q1_right_gap;
                } else {
                    p2f += score_model.get_gap_penalty(q1_right_gap + q2_right_gap)
                        + score_model.inv_ext_p * q1_right_gap;
                }

                qcons2f += q1_right_gap + q2_right_gap;
            }
            p2f += (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * get_q2_fw(query_2).size();

            auto [p2b, cigar2b, rcons2b, qcons2b] = align_segment(
                aligner, score_model, get_q2_bw(query_rc_2), get_t2_bw(target_2),
                penalty_to_score,
                heuristics_length_cutoff
            );

            Offset r_consumed = rcons2f + rcons2b;
            assert(r_consumed == target_2.size());

            Offset q_consumed = qcons2f + qcons2b;
            assert(q_consumed == query_2.size() + query_rc_2.size());

            return std::make_tuple(
                p2f + p2b,
                std::move(cigar2f + std::move(cigar2b)),
                r_consumed,
                q_consumed
            );
        }
    }

    bool all_on_left(Offset t1_max, Offset q1_max, Offset t2_max, Offset q2_max) const {
        return q1 + q1_left_gap == q1_max
                            && q2 + q2_right_gap == 0
                            && t2 + t2_left_gap == t2_max && (t1_max > 0 || q1_max > 0 || q2_max > 0);
    }

    std::tuple<Score, std::string, Offset, Offset>
    align_all_left(wfa::WFAligner& aligner,
                   const ScoreModel &score_model,
                   std::string_view target_1,
                   std::string_view query_1,
                   std::string_view query_rc_1,
                   bool get_cigar,
                   SOffset heuristics_length_cutoff,
                   bool penalty_to_score) const {
        if (!t1_left_gap && !t1_right_gap && !q1_left_gap && !q2_left_gap) {
            std::string query_1_cat(query_1);
            query_1_cat += query_rc_1;

            auto [p, cigar, r_consumed, q_consumed] = align_segment(
                aligner, score_model, query_1_cat, target_1,
                penalty_to_score,
                heuristics_length_cutoff
            );

            return std::make_tuple(
                p + (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * query_rc_1.size(),
                std::move(cigar),
                r_consumed,
                q_consumed
            );
        } else {
            auto [p1f, cigar1f, rcons1f, qcons1f] = align_segment(
                aligner, score_model, get_q1_fw(query_1), get_t1_fw(target_1),
                penalty_to_score,
                heuristics_length_cutoff
            );
            if (t1_left_gap + t1_right_gap) {
                cigar1f += std::to_string(t1_left_gap + t1_right_gap) + "I";
                p1f += penalty_to_score
                    ? score_model.get_gap_score(t1_left_gap + t1_right_gap)
                    : score_model.get_gap_penalty(t1_left_gap + t1_right_gap);
                rcons1f += t1_left_gap + t1_right_gap;
            }
            if (q1_left_gap + q2_left_gap) {
                cigar1f += std::to_string(q1_left_gap + q2_left_gap) + "D";
                if (penalty_to_score) {
                    p1f += score_model.get_gap_score(q1_left_gap + q2_left_gap)
                        + score_model.inv_ext_s * q2_left_gap;
                } else {
                    p1f += score_model.get_gap_penalty(q1_left_gap + q2_left_gap)
                        + score_model.inv_ext_p * q2_left_gap;
                }

                qcons1f += q1_left_gap + q2_left_gap;
            }
            auto [p1b, cigar1b, rcons1b, qcons1b] = align_segment(
                aligner, score_model, get_q1_bw(query_rc_1), get_t1_bw(target_1),
                penalty_to_score,
                heuristics_length_cutoff
            );
            p1b += (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * get_q1_bw(query_rc_1).size();

            Offset r_consumed = rcons1f + rcons1b;
            assert(r_consumed == target_1.size());

            Offset q_consumed = qcons1f + qcons1b;
            assert(q_consumed == query_1.size() + query_rc_1.size());

            return std::make_tuple(
                p1f + p1b,
                std::move(cigar1f + std::move(cigar1b)),
                r_consumed,
                q_consumed
            );
        }
    }
};

inline std::ostream& operator<<(std::ostream& out, const Breakpoints &bp) {
    out << bp.min_p << "\t"
        << "t1: " << bp.t1 << " " << bp.t1_left_gap << "," << bp.t1_right_gap << "\t"
        << "q1: " << bp.q1 << " " << bp.q1_left_gap << "," << bp.q1_right_gap << "\t"
        << "t2: " << bp.t2 << " " << bp.t2_left_gap << "," << bp.t2_right_gap << "\t"
        << "q2: " << bp.q2 << " " << bp.q2_left_gap << "," << bp.q2_right_gap << "\n";
    return out;
}

template <bool inv_open, bool inv_ext, bool is_fwd, bool is_open = (is_fwd != inv_ext)>
void update_breakpoints(wfa::WFAligner& aligner,
                        const ScoreModel &score_model,
                        Breakpoints &breakpoints,
                        WFAIterator<false, false, true> &wfa_it_1_fwd_in,
                        WFAIterator<true, true, false> &wfa_it_1_bwd_in,
                        WFAIterator<false, true, true> &wfa_it_2_fwd_in,
                        WFAIterator<false, false, false> &wfa_it_2_bwd_in,
                        Offset q,
                        Offset t,
                        Offset ext,
                        bool exhausted,
                        SOffset heuristics_length_cutoff) {
    assert(wfa_it_1_fwd_in.max_q() == wfa_it_2_fwd_in.max_q());
    assert(wfa_it_1_bwd_in.max_q() == wfa_it_2_bwd_in.max_q());
    assert(wfa_it_1_fwd_in.max_t() == wfa_it_1_bwd_in.max_t());
    assert(wfa_it_2_fwd_in.max_t() == wfa_it_2_bwd_in.max_t());

    const SOffset q1_max = wfa_it_1_fwd_in.max_q();
    const SOffset q2_max = wfa_it_1_bwd_in.max_q();
    const SOffset t1_max = wfa_it_1_fwd_in.max_t();
    const SOffset t2_max = wfa_it_2_fwd_in.max_t();

    auto update = [&aligner,&breakpoints,&score_model,q1_max,q2_max,t1_max,t2_max,
                   q,t,ext,
                   &wfa_it_1_fwd_in,&wfa_it_1_bwd_in,
                   &wfa_it_2_fwd_in,&wfa_it_2_bwd_in,&exhausted,
                   &heuristics_length_cutoff](SOffset q1, SOffset q1p, SOffset q2, SOffset q2p,
                                                             SOffset t1, SOffset t1p, SOffset t2, SOffset t2p,
                                                             Penalty update_p,
                                                             SOffset q1_left_gap = 0,
                                                             SOffset q1_right_gap = 0,
                                                             SOffset q2_left_gap = 0,
                                                             SOffset q2_right_gap = 0,
                                                             SOffset t1_left_gap = 0,
                                                             SOffset t1_right_gap = 0,
                                                             SOffset t2_left_gap = 0,
                                                             SOffset t2_right_gap = 0) {
        assert(q1 >= 0);
        assert(q1p >= 0);
        assert(q2 >= 0);
        assert(q2p >= 0);
        assert(t1 >= 0);
        assert(t1p >= 0);
        assert(t2 >= 0);
        assert(t2p >= 0);
        assert(q1_left_gap >= 0);
        assert(q2_left_gap >= 0);
        assert(t1_left_gap >= 0);
        assert(t2_left_gap >= 0);
        assert(q1_right_gap >= 0);
        assert(q2_right_gap >= 0);
        assert(t1_right_gap >= 0);
        assert(t2_right_gap >= 0);
        assert(t1 + t1p + t1_left_gap + t1_right_gap == t1_max);
        assert(t2 + t2p + t2_left_gap + t2_right_gap  == t2_max);
        assert(q1 + q1p + q1_left_gap + q1_right_gap  == q1_max);
        assert(q2 + q2p + q2_left_gap + q2_right_gap  == q2_max);

        if (update_p >= breakpoints.min_p)
            return;

        Breakpoints bp_new;
        bp_new.t1 = t1;
        bp_new.q1 = q1;
        bp_new.q2 = q2;
        bp_new.t2 = t2;
        bp_new.t1_left_gap = t1_left_gap;
        bp_new.t1_right_gap = t1_right_gap;
        bp_new.t2_left_gap = t2_left_gap;
        bp_new.t2_right_gap = t2_right_gap;
        bp_new.q1_left_gap = q1_left_gap;
        bp_new.q1_right_gap = q1_right_gap;
        bp_new.q2_left_gap = q2_left_gap;
        bp_new.q2_right_gap = q2_right_gap;

        if (bp_new.all_on_right(t1_max, q1_max, t2_max, t2_max)) {
            // t1 will align to nothing, t2 will align to the concatenation of query_2 and query_rc_2
            if (t2_max == 0) {
                // penalty already covers query_rc_2, now we need to cover query_2
                update_p += score_model.get_gap_penalty(q1_max)
                            + score_model.inv_ext_p * q1_max;
            } else {
                auto [newp, dummy_cigar, rcons, qcons] = bp_new.align_all_right(
                    aligner,
                    score_model,
                    wfa_it_2_fwd_in.get_target(),
                    wfa_it_2_fwd_in.get_query(),
                    wfa_it_2_bwd_in.get_query(),
                    false,
                    heuristics_length_cutoff,
                    false
                );
                update_p = newp;
            }

            if (update_p >= breakpoints.min_p)
                return;
        } else if (bp_new.all_on_left(t1_max, q1_max, t2_max, t2_max)) {
            // t1 will align to the concatenation of query_1 and query_rc_1, t2 will align to nothing
            if (t1_max == 0) {
                // penalty already covers query_1, now we need to cover query_rc_1
                update_p += score_model.get_gap_penalty(q2_max)
                            + score_model.inv_ext_p * q2_max;
            } else {
                auto [newp, dummy_cigar, rcons, qcons] = bp_new.align_all_left(
                    aligner,
                    score_model,
                    wfa_it_1_fwd_in.get_target(),
                    wfa_it_1_fwd_in.get_query(),
                    wfa_it_1_bwd_in.get_query(),
                    false,
                    heuristics_length_cutoff,
                    false
                );
                update_p = newp;
            }

            if (update_p >= breakpoints.min_p)
                return;
        }

        bp_new.min_p = update_p;
        breakpoints = bp_new;

        std::cerr << "FOO\t" << breakpoints << "\n";
    };

    auto &wfa_it_1_fwd = [&]() -> auto& {
        if constexpr (is_open && is_fwd) {
            return wfa_it_1_fwd_in;
        } else if constexpr (is_open) {
            // swap fw and bw
            return wfa_it_1_bwd_in;
        } else if constexpr (is_fwd) {
            // swap wfa_1 and wfa_2
            return wfa_it_2_fwd_in;
        } else {
            // swap wfa_1 and wfa_2, and fw and bw
            return wfa_it_2_bwd_in;
        }
    }();

    auto &wfa_it_1_bwd = [&]() -> auto& {
        if constexpr (is_open && is_fwd) {
            return wfa_it_1_bwd_in;
        } else if constexpr (is_open) {
            // swap fw and bw
            return wfa_it_1_fwd_in;
        } else if constexpr (is_fwd) {
            // swap wfa_1 and wfa_2
            return wfa_it_2_bwd_in;
        } else {
            // swap wfa_1 and wfa_2, and fw and bw
            return wfa_it_2_fwd_in;
        }
    }();

    auto &wfa_it_2_fwd = [&]() -> auto& {
        if constexpr (is_open && is_fwd) {
            return wfa_it_2_fwd_in;
        } else if constexpr (is_open) {
            // swap fw and bw
            return wfa_it_2_bwd_in;
        } else if constexpr (is_fwd) {
            // swap wfa_1 and wfa_2
            return wfa_it_1_fwd_in;
        } else {
            // swap wfa_1 and wfa_2, and fw and bw
            return wfa_it_1_bwd_in;
        }
    }();

    auto &wfa_it_2_bwd = [&]() -> auto& {
        if constexpr (is_open && is_fwd) {
            return wfa_it_2_bwd_in;
        } else if constexpr (is_open) {
            // swap fw and bw
            return wfa_it_2_fwd_in;
        } else if constexpr (is_fwd) {
            // swap wfa_1 and wfa_2
            return wfa_it_1_bwd_in;
        } else {
            // swap wfa_1 and wfa_2, and fw and bw
            return wfa_it_1_fwd_in;
        }
    }();

    assert(wfa_it_1_fwd.max_t() == wfa_it_1_bwd.max_t());
    assert(wfa_it_2_fwd.max_t() == wfa_it_2_bwd.max_t());
    assert(wfa_it_1_fwd.max_q() == wfa_it_2_fwd.max_q());
    assert(wfa_it_1_bwd.max_q() == wfa_it_2_bwd.max_q());

    auto permute_then_update = [&update,&wfa_it_1_fwd,&wfa_it_1_bwd,&wfa_it_2_fwd,&wfa_it_2_bwd](SOffset q1, SOffset q1p, SOffset q2, SOffset q2p,
                                         SOffset t1, SOffset t1p, SOffset t2, SOffset t2p,
                                         Penalty update_p,
                                         SOffset q1_left_gap = 0,
                                        SOffset q1_right_gap = 0,
                                        SOffset q2_left_gap = 0,
                                        SOffset q2_right_gap = 0,
                                        SOffset t1_left_gap = 0,
                                        SOffset t1_right_gap = 0,
                                        SOffset t2_left_gap = 0,
                                        SOffset t2_right_gap = 0) {
        assert(t1 + t1p + t1_left_gap + t1_right_gap == wfa_it_1_fwd.max_t());
        assert(t2 + t2p + t2_left_gap + t2_right_gap == wfa_it_2_fwd.max_t());
        assert(q1 + q1p + q1_left_gap + q1_right_gap == wfa_it_1_fwd.max_q());
        assert(q2 + q2p + q2_left_gap + q2_right_gap == wfa_it_1_bwd.max_q());

        if constexpr (is_open && is_fwd) {
            update(q1, q1p, q2, q2p, t1, t1p, t2, t2p, update_p,
                   q1_left_gap, q1_right_gap,
                   q2_left_gap, q2_right_gap,
                   t1_left_gap, t1_right_gap,
                   t2_left_gap, t2_right_gap);
        } else if constexpr (is_open && !is_fwd) {
            // we calculated t1' as t1 and t2' as t2
            // we calculated q2' as q1 and q1' as q2
            update(q2p, q2, q1p, q1, t1p, t1, t2p, t2, update_p,
                   q2_right_gap, q2_left_gap,
                   q1_right_gap, q1_left_gap,
                   t1_right_gap, t1_left_gap,
                   t2_right_gap, t2_left_gap);
        } else if constexpr (!is_open && is_fwd) {
            // we calculated t2 as t1 and t1 as t2
            // we calculate q1' as q1 and q2' as q2
            update(q1p, q1, q2p, q2, t2, t2p, t1, t1p, update_p,
                   q1_right_gap, q1_left_gap,
                   q2_right_gap, q2_left_gap,
                   t2_left_gap, t2_right_gap,
                   t1_left_gap, t1_right_gap);
        } else {
            // we calculated t2' as t1 and t2' as t2
            // we calculate q2 as q1 and q1 as q2
            update(q2, q2p, q1, q1p, t2p, t2, t1p, t1, update_p,
                   q2_left_gap, q2_right_gap,
                   q1_left_gap, q1_right_gap,
                   t2_right_gap, t2_left_gap,
                   t1_right_gap, t1_left_gap);
        }
    };

    SOffset q1_b = q;
    SOffset q1_e = q + ext + 1;

    SOffset t1_b = t;
    SOffset t1_e = t + ext + 1;

    // if (exhausted) {
    //     assert(q1_e - 1 == wfa_it_1_fwd.max_q() || t1_e - 1 == wfa_it_1_fwd.max_t());
    //     if (q1_e - 1 == wfa_it_1_fwd.max_q()) {
    //         SOffset q1 = q1_e - 1;
    //         SOffset q1p = 0;
    //         SOffset t1 = t1_e - 1;
    //     }
    //     // SOffset q1 = q1_e - 1;
    //     // SOffset q1p = wfa_it_1_fwd.max_q() - q1;
    //     // SOffset t1 = t1_e - 1;
    //     // SOffset t1p = wfa_it_1_fwd.max_t() - t1;

    //     // SOffset t2p = std::get<1>(wfa_it_2_bwd.get_frp());
    //     // SOffset q2 = std::get<0>(wfa_it_2_bwd.get_frp());
    //     // SOffset t2 = wfa_it_2_fwd.max_t() - t2p;
    //     // SOffset q2p = wfa_it_2_bwd.max_q() - q2;

    //     // permute_then_update(q1, q1p, q2, q2p, t1, t1p, t2, t2p, 0);
    //     return;
    // }

    auto [q1, t1, p1, ext_1f] = wfa_it_1_fwd_in.get_frp();
    auto [q2p, t1p, p1p, ext_1b] = wfa_it_1_bwd_in.get_frp();
    auto [q1p, t2, p2, ext_2f] = wfa_it_2_fwd_in.get_frp();
    auto [q2, t2p, p2p, ext_2b] = wfa_it_2_bwd_in.get_frp();

    SOffset q1_gap = q1_max - q1 - q1p;
    SOffset q2_gap = q2_max - q2 - q2p;
    SOffset t1_gap = t1_max - t1 - t1p;
    SOffset t2_gap = t2_max - t2 - t2p;

    Penalty cur_p;
    if constexpr (!inv_open && !inv_ext && is_fwd) {
        cur_p = wfa_it_1_fwd_in.get_p();
    } else if constexpr (inv_open && inv_ext && !is_fwd) {
        cur_p = wfa_it_1_bwd_in.get_p();
    } else if constexpr (!inv_open && inv_ext && is_fwd) {
        cur_p = wfa_it_2_fwd_in.get_p();
    } else {
        cur_p = wfa_it_2_bwd_in.get_p();
    }

    if (q1_gap > 0 || q2_gap > 0 || t1_gap > 0 || t2_gap > 0) {
        SOffset ext_q_1 = ext_1f + ext_2f;
        SOffset ext_q_2 = ext_1b + ext_2b;
        SOffset ext_t_1 = ext_1f + ext_1b;
        SOffset ext_t_2 = ext_2f + ext_2b;

        if (-q1_gap <= ext_q_1 && -q2_gap <= ext_q_2 && -t1_gap <= ext_t_1 && -t2_gap <= ext_t_2) {
            assert(q1 + q1p + q1_gap == q1_max);
            assert(q2 + q2p + q2_gap == q2_max);
            assert(t1 + t1p + t1_gap == t1_max);
            assert(t2 + t2p + t2_gap == t2_max);

            if (q1_gap < 0) {
                q1p += q1_gap;
                t1p += q1_gap;
                t1_gap -= q1_gap;

                q1_gap = 0;

                if (t1p < 0 || t1_gap < 0 || q1p < 0 || q2_gap < 0)
                    return;

                assert(q1p >= 0);
                assert(q1 + q1p == q1_max);
                assert(t1 + t1p + t1_gap == t1_max);
            }

            if (q2_gap < 0) {
                q2p += q2_gap;
                t2 += q2_gap;
                t2_gap -= q2_gap;

                q2_gap = 0;

                if (t2 < 0 || t2_gap < 0 || q2p < 0 || q1_gap < 0)
                    return;

                assert(q2 + q2p == q2_max);
                assert(t2 + t2p + t2_gap == t2_max);
            }

            if (t1_gap < 0 || t2_gap < 0)
                return;

            assert(q1_gap >= 0);
            assert(q2_gap >= 0);

            Penalty p = p1 + p1p + p2 + p2p
                            + score_model.get_gap_penalty(q1_gap)
                            + score_model.get_gap_penalty(q2_gap)
                            + score_model.get_gap_penalty(t1_gap)
                            + score_model.get_gap_penalty(t2_gap);

            SOffset t1_left_gap = t1_gap;
            SOffset t1_right_gap = 0;
            SOffset t2_left_gap = 0;
            SOffset t2_right_gap = t2_gap;
            SOffset q1_left_gap = q1_gap;
            SOffset q1_right_gap = 0;
            SOffset q2_left_gap = 0;
            SOffset q2_right_gap = q2_gap;

            update(q1, q1p,
                   q2, q2p,
                   t1, t1p,
                   t2, t2p,
                   p,
                   q1_left_gap, q1_right_gap,
                   q2_left_gap, q2_right_gap,
                   t1_left_gap, t1_right_gap,
                   t2_left_gap, t2_right_gap);
        }

        return;
    }

    // we have t1 and q1 (and diag_1_f)
    // we then compute a t2 (and diag_2_f) that works using q1'
    // we then compute a q2 (and diag_2_b) that works using t2'

    // alternately
    // we then compute a q2' (and diag_1_b) that works using t1'
    // we then compute a t2' (and diag_2_b) that works using q2

    Diag diag_1_f = WFAIterator<>::get_diag(q1_b, t1_b);
    assert(wfa_it_1_fwd.get_min_diag() <= diag_1_f);
    assert(diag_1_f <= wfa_it_1_fwd.get_max_diag());

    if (wfa_it_1_fwd.get_global_min(diag_1_f) >= breakpoints.min_p) {
        wfa_it_1_fwd.disable_diag(diag_1_f);
        return;
    }

    auto [it_b, it_e] = wfa_it_1_fwd.get_min(diag_1_f, q1_b, q1_e);
    assert(it_b != it_e);

    Penalty min_update_p_1 = it_b->penalty;
    if (min_update_p_1 >= breakpoints.min_p)
        return;

    SOffset q1p_b = wfa_it_1_fwd.max_q() - (q1_e - 1);
    SOffset q1p_e = wfa_it_1_fwd.max_q() - q1_b + 1;
    q1p_e = std::min(q1p_e, wfa_it_2_fwd.get_global_end() + 1);
    if (q1p_b >= q1p_e)
        return;

    SOffset t1p_b = wfa_it_1_fwd.max_t() - (t1_e - 1);
    SOffset t1p_e = wfa_it_1_fwd.max_t() - t1_b + 1;
    t1p_e = std::min(t1p_e, wfa_it_1_bwd.get_global_r_end() + 1);
    if (t1p_b >= t1p_e)
        return;

    // now, try to find a diag in wfa_it_2_fwd that overlaps q1' so we can compute q1, t1, and t2
    auto d_2f = wfa_it_2_fwd.get_overlaps_with_q(q1p_b, q1p_e);
    if (d_2f.empty())
        return;

    auto d_1b = wfa_it_1_bwd.get_overlaps_with_t(t1p_b, t1p_e);

    // std::cerr << "foocomp\t" << dir << "\t" << d_2f.size() << " vs. " << d_1b.size() << "\n";

    if (d_2f.size() <= d_1b.size()) {
        for (auto [q1p, t2] : d_2f) {
            Diag diag_2_f = wfa_it_2_fwd.get_diag(q1p, t2);
            Penalty min_update_p_2 = wfa_it_2_fwd.get_global_min(diag_2_f);
            if (min_update_p_2 >= breakpoints.min_p) {
                wfa_it_2_fwd.disable_diag(diag_2_f);
                continue;
            }

            min_update_p_2 += min_update_p_1;
            if (min_update_p_2 >= breakpoints.min_p)
                continue;

            auto [kt_b, kt_e] = wfa_it_2_fwd.get_min(diag_2_f, q1p, q1p + 1);
            if (kt_b == kt_e)
                continue;

            min_update_p_2 = min_update_p_1 + kt_b->penalty;
            if (min_update_p_2 >= breakpoints.min_p)
                continue;

            SOffset t2p = wfa_it_2_fwd.max_t() - t2;

            SOffset q1 = wfa_it_1_fwd.max_q() - q1p;
            assert(q1_b <= q1);
            assert(q1 < q1_e);

            SOffset t1 = diag_1_f + q1;
            assert(t1_b <= t1);
            assert(t1 < t1_e);

            auto [it_b, it_e] = wfa_it_1_fwd.get_min(diag_1_f, q1, q1 + 1);
            if (it_b == it_e)
                break;

            // we now have q1 and t1
            min_update_p_2 = it_b->penalty + kt_b->penalty;
            if (min_update_p_2 >= breakpoints.min_p)
                continue;

            // we now have q1, t1, and t2

            // diag_2_b = t2' - q2'' = t2' - q2 => q2 = t2' - diag_2_b
            auto d_2b = wfa_it_2_bwd.get_overlaps_with_t(t2p, t2p + 1);
            if (d_2b.empty())
                continue;

            SOffset q2 = d_2b.front().first;
            Diag diag_2_b = wfa_it_2_bwd.get_diag(q2, t2p);
            Penalty min_update_p_3 = wfa_it_2_bwd.get_global_min(diag_2_b);
            if (min_update_p_3 >= breakpoints.min_p) {
                wfa_it_2_bwd.disable_diag(diag_2_b);
                continue;
            }

            min_update_p_3 += min_update_p_2;
            if (min_update_p_3 >= breakpoints.min_p)
                continue;

            auto [lt_b, lt_e] = wfa_it_2_bwd.get_min(diag_2_b, q2, q2 + 1);
            if (lt_b == lt_e)
                continue;

            min_update_p_3 = min_update_p_2 + lt_b->penalty;
            if (min_update_p_3 >= breakpoints.min_p)
                continue;

            SOffset q2p = wfa_it_2_bwd.max_q() - q2;
            SOffset t1p = wfa_it_1_fwd.max_t() - t1;

            // we now have q1, t1, t2, q2!
            // find the corresponding jt_b, jt_e
            Diag diag_1_b = t1p - q2p;
            if (wfa_it_1_bwd.get_global_min(diag_1_b) >= breakpoints.min_p) {
                wfa_it_1_bwd.disable_diag(diag_1_b);
                continue;
            }

            auto [jt_b, jt_e] = wfa_it_1_bwd.get_min(diag_1_b, q2p, q2p + 1);
            if (jt_b == jt_e)
                continue;

            Penalty update_p = min_update_p_3 + jt_b->penalty;
            permute_then_update(q1, q1p, q2, q2p, t1, t1p, t2, t2p, update_p);
        }
    } else {
        for (auto [q2p, t1p] : d_1b) {
            Diag diag_1_b = wfa_it_1_bwd.get_diag(q2p, t1p);
            Penalty min_update_p_2 = wfa_it_1_bwd.get_global_min(diag_1_b);
            if (min_update_p_2 >= breakpoints.min_p) {
                wfa_it_1_bwd.disable_diag(diag_1_b);
                continue;
            }

            min_update_p_2 += min_update_p_1;
            if (min_update_p_2 >= breakpoints.min_p)
                continue;

            auto [jt_b, jt_e] = wfa_it_1_bwd.get_min(diag_1_b, q2p, q2p + 1);
            if (jt_b == jt_e)
                continue;

            min_update_p_2 = min_update_p_1 + jt_b->penalty;
            if (min_update_p_2 >= breakpoints.min_p)
                continue;

            // we now have t1, q1, and q2
            SOffset t1 = wfa_it_1_fwd.max_t() - t1p;
            SOffset q2 = wfa_it_1_bwd.max_q() - q2p;

            SOffset q1 = t1 - diag_1_f;
            auto [it_b, it_e] = wfa_it_1_fwd.get_min(diag_1_f, q1, q1 + 1);
            if (it_b == it_e)
                continue;

            SOffset q1p = wfa_it_1_fwd.max_q() - q1;

            // now find t2'
            // diag_2_b = t2' - q2
            auto d_2b = wfa_it_2_bwd.get_overlaps_with_q(q2, q2 + 1);
            if (d_2b.empty())
                continue;

            SOffset t2p = d_2b.front().second;
            Diag diag_2_b = wfa_it_2_bwd.get_diag(q2, t2p);
            Penalty min_update_p_3 = wfa_it_2_bwd.get_global_min(diag_2_b);
            if (min_update_p_3 >= breakpoints.min_p) {
                wfa_it_2_bwd.disable_diag(diag_2_b);
                continue;
            }

            min_update_p_3 += min_update_p_2;
            if (min_update_p_2 >= breakpoints.min_p)
                continue;

            auto [lt_b, lt_e] = wfa_it_2_bwd.get_min(diag_2_b, q2, q2 + 1);
            if (lt_b == lt_e)
                continue;

            // we now have t1, q1, q2 and t2
            SOffset t2 = wfa_it_2_bwd.max_t() - t2p;

            Diag diag_2_f = wfa_it_2_fwd.get_diag(q1p, t2);

            auto [kt_b, kt_e] = wfa_it_2_fwd.get_min(diag_2_f, q1p, q1p + 1);
            if (kt_b == kt_e)
                continue;

            Penalty update_p = it_b->penalty + jt_b->penalty + kt_b->penalty + lt_b->penalty;
            permute_then_update(q1, q1p, q2, q2p, t1, t1p, t2, t2p, update_p);
        }
    }
}

std::tuple<Score, std::string, SOffset, SOffset, SOffset, SOffset,
           Score, std::string, SOffset, SOffset, SOffset, SOffset>
get_alignment_cigars(wfa::WFAligner& aligner,
                     const ScoreModel &score_model,
                     std::string_view query_1,
                     std::string_view query_rc_1,
                     std::string_view target_1,
                     std::string_view query_2,
                     std::string_view query_rc_2,
                     std::string_view target_2,
                     const Breakpoints &breakpoints,
                     SOffset heuristics_length_cutoff) {
    assert(breakpoints.min_p < ScoreModel::inf_p);
    assert(breakpoints.q1 <= query_1.size());
    assert(breakpoints.q2 <= query_rc_1.size());
    assert(breakpoints.t1 <= target_1.size());
    assert(breakpoints.t2 <= target_2.size());

    auto align_segment = [&aligner,&score_model,&heuristics_length_cutoff](std::string_view query,
                                                 std::string_view target,
                                                 Diag min_k = min_diag,
                                                 Diag max_k = max_diag) -> std::tuple<Score, std::string, Offset, Offset> {
        return Breakpoints::align_segment(
            aligner, score_model,
            query, target,
            true,
            heuristics_length_cutoff,
            min_diag, max_diag
        );
    };

    Score score_1_fw = 0;
    Score score_1_bw = 0;
    Score score_2_fw = 0;
    Score score_2_bw = 0;
    Offset r_consumed_1 = 0;
    Offset r_consumed_1_bw = 0;
    Offset q_consumed_1 = 0;
    Offset q_consumed_1_bw = 0;
    Offset r_consumed_2 = 0;
    Offset r_consumed_2_bw = 0;
    Offset q_consumed_2 = 0;
    Offset q_consumed_2_bw = 0;

    SOffset inv_length_1 = 0;
    SOffset inv_length_r_1 = 0;
    SOffset inv_length_2 = 0;
    SOffset inv_length_r_2 = 0;

    std::string cigar_1;
    std::string cigar_1_bw;
    std::string cigar_2;
    std::string cigar_2_bw;

    std::string query_1_cat;
    std::string query_2_cat;

    Score score_1 = score_model.inv_open_s;
    Score score_2 = 0;

    if (breakpoints.all_on_right(target_1.size(), query_1.size(), target_2.size(), query_rc_1.size())) {
        // t1 will align to nothing, t2 will align to the concatenation of query_2 and query_rc_2
        if (target_1.size()) {
            cigar_1 = std::to_string(target_1.size()) + "I";
            r_consumed_1 = target_1.size();
            score_1_fw = score_model.get_gap_score(target_1.size());
        }

        std::tie(score_2_fw, cigar_2, r_consumed_2, q_consumed_2) = breakpoints.align_all_right(
            aligner,
            score_model,
            target_2,
            query_2,
            query_rc_2,
            true,
            heuristics_length_cutoff,
            true
        );
        inv_length_2 = query_2.size();
        inv_length_r_2 = cigar_get_target_pos(cigar_2, query_2.size());
    } else if (breakpoints.all_on_left(target_1.size(), query_1.size(), target_2.size(), query_rc_1.size())) {
        if (target_2.size()) {
            cigar_2 = std::to_string(target_2.size()) + "I";
            r_consumed_2 = target_2.size();
            score_2_bw = score_model.get_gap_score(target_2.size());
        }
        std::tie(score_1_fw, cigar_1, r_consumed_1, q_consumed_1) = breakpoints.align_all_left(
            aligner,
            score_model,
            target_1,
            query_1,
            query_rc_1,
            true,
            heuristics_length_cutoff,
            true
        );
        inv_length_1 = query_rc_1.size();
        inv_length_r_1 = target_1.size() - cigar_get_target_pos(cigar_1, query_1.size());
    } else {
        std::tie(score_1_fw, cigar_1, r_consumed_1, q_consumed_1) = align_segment(
            breakpoints.get_q1_fw(query_1),
            breakpoints.get_t1_fw(target_1)
            // wfa_it_1_fwd.get_min_diag(),
            // wfa_it_1_fwd.get_max_diag()
        );

        if (breakpoints.t1_left_gap > 0) {
            cigar_1 += std::to_string(breakpoints.t1_left_gap) + "I";
            score_1_fw += score_model.get_gap_score(breakpoints.t1_left_gap);
            r_consumed_1 += breakpoints.t1_left_gap;
        }

        if (breakpoints.q1_left_gap > 0) {
            cigar_1 += std::to_string(breakpoints.q1_left_gap) + "D";
            score_1_fw += score_model.get_gap_score(breakpoints.q1_left_gap);
            q_consumed_1 += breakpoints.q1_left_gap;
        }

        std::tie(score_1_bw, cigar_1_bw, r_consumed_1_bw, q_consumed_1_bw) = align_segment(
            breakpoints.get_q1_bw(query_rc_1),
            breakpoints.get_t1_bw(target_1)
        );

        inv_length_1 = breakpoints.get_q1_bw(query_rc_1).size();
        inv_length_r_1 = breakpoints.get_t1_bw(target_1).size();

        if (breakpoints.t1_right_gap > 0) {
            cigar_1_bw = std::to_string(breakpoints.t1_right_gap) + "I" + cigar_1_bw;
            score_1_bw += score_model.get_gap_score(breakpoints.t1_right_gap);
            inv_length_r_1 += breakpoints.t1_right_gap;
            r_consumed_1_bw += breakpoints.t1_right_gap;
        }

        if (breakpoints.q2_left_gap > 0) {
            cigar_1_bw = std::to_string(breakpoints.q2_left_gap) + "D" + cigar_1_bw;
            score_1_bw += score_model.get_gap_score(breakpoints.q2_left_gap);
            inv_length_1 += breakpoints.q2_left_gap;
            q_consumed_1_bw += breakpoints.q2_left_gap;
        }

        std::tie(score_2_fw, cigar_2, r_consumed_2, q_consumed_2) = align_segment(
            breakpoints.get_q2_fw(query_2),
            breakpoints.get_t2_fw(target_2)
            // wfa_it_2_fwd.get_min_diag(),
            // wfa_it_2_fwd.get_max_diag()
        );

        inv_length_2 = breakpoints.get_q2_fw(query_2).size();
        inv_length_r_2 = breakpoints.get_t2_fw(target_2).size();

        if (breakpoints.t2_left_gap > 0) {
            cigar_2 += std::to_string(breakpoints.t2_left_gap) + "I";
            score_2_fw += score_model.get_gap_score(breakpoints.t2_left_gap);
            inv_length_r_2 += breakpoints.t2_left_gap;
            r_consumed_2 += breakpoints.t2_left_gap;
        }

        if (breakpoints.q1_right_gap > 0) {
            cigar_2 += std::to_string(breakpoints.q1_right_gap) + "D";
            score_2_fw += score_model.get_gap_score(breakpoints.q1_right_gap);
            inv_length_2 += breakpoints.q1_right_gap;
            q_consumed_2 += breakpoints.q1_right_gap;
        }

        std::tie(score_2_bw, cigar_2_bw, r_consumed_2_bw, q_consumed_2_bw) = align_segment(
            breakpoints.get_q2_bw(query_rc_2),
            breakpoints.get_t2_bw(target_2)
        );

        if (breakpoints.t2_right_gap > 0) {
            cigar_2_bw = std::to_string(breakpoints.t2_right_gap) + "I" + cigar_2_bw;
            score_2_bw += score_model.get_gap_score(breakpoints.t2_right_gap);
            r_consumed_2_bw += breakpoints.t2_right_gap;
        }

        if (breakpoints.q2_right_gap > 0) {
            cigar_2_bw = std::to_string(breakpoints.q2_right_gap) + "D" + cigar_2_bw;
            score_2_bw += score_model.get_gap_score(breakpoints.q2_right_gap);
            q_consumed_2_bw += breakpoints.q2_right_gap;
        }

        score_1 += inv_length_1 * score_model.inv_ext_s;
        score_2 += inv_length_2 * score_model.inv_ext_s;
    }

    r_consumed_1 += r_consumed_1_bw;
    q_consumed_1 += q_consumed_1_bw;

    r_consumed_2 += r_consumed_2_bw;
    q_consumed_2 += q_consumed_2_bw;

    assert(r_consumed_1 + r_consumed_2 == target_1.size() + target_2.size());
    assert(q_consumed_1 + q_consumed_2 == query_1.size() + query_rc_1.size());

    cigar_1 += std::move(cigar_1_bw);
    cigar_2 += std::move(cigar_2_bw);

    score_1 += score_1_fw + score_1_bw;
    score_2 += score_2_fw + score_2_bw;

    return std::make_tuple(score_1, std::move(cigar_1), r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
                           score_2, std::move(cigar_2), r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2);
}

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
              SOffset heuristics_length_cutoff,
              SOffset min_wavefront_length,
              SOffset max_distance_threshold) {
    assert(query_1.size() == query_2.size());
    assert(query_rc_1.size() == query_rc_2.size());
    assert(query_1 == reverse_complement(query_2));
    assert(query_rc_1 == reverse_complement(query_rc_2));

    Breakpoints breakpoints {
        .min_p = ScoreModel::inf_p,
        .t1 = target_1.size() + 1,
        .t2 = target_2.size() + 1,
        .q1 = query_1.size() + 1,
        .q2 = query_rc_1.size() + 1
    };

    const bool use_heuristics = ((min_wavefront_length < max_diag_width)
                                || (max_distance_threshold < max_offset)) && false;

    WFAIterator<false, false, true> wfa_it_1_fwd(
        score_model,
        query_1.begin(), query_1.end(),
        target_1.begin(), target_1.end(),
        min_wavefront_length,
        max_distance_threshold
    );

    WFAIterator<true, true, false> wfa_it_1_bwd(
        score_model,
        query_rc_1.begin(), query_rc_1.end(),
        target_1.begin(), target_1.end(),
        min_wavefront_length,
        max_distance_threshold
    );

    WFAIterator<false, true, true> wfa_it_2_fwd(
        score_model,
        query_2.begin(), query_2.end(),
        target_2.begin(), target_2.end(),
        min_wavefront_length,
        max_distance_threshold
    );

    WFAIterator<false, false, false> wfa_it_2_bwd(
        score_model,
        query_rc_2.begin(), query_rc_2.end(),
        target_2.begin(), target_2.end(),
        min_wavefront_length,
        max_distance_threshold
    );

    assert(wfa_it_1_fwd.get_p() == wfa_it_1_bwd.get_p());
    assert(wfa_it_1_bwd.get_p() == wfa_it_2_fwd.get_p());
    assert(wfa_it_2_fwd.get_p() == wfa_it_2_bwd.get_p());

    update_breakpoints<false, false, true>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, 0, 0, 0, false, heuristics_length_cutoff);

    bool exhausted_1_fwd = target_1.empty() || query_1.empty();
    bool exhausted_1_bwd = target_1.empty() || query_rc_1.empty();
    bool exhausted_2_fwd = target_2.empty() || query_2.empty();
    bool exhausted_2_bwd = target_2.empty() || query_rc_2.empty();

    while (!exhausted_1_fwd || !exhausted_1_bwd || !exhausted_2_fwd || !exhausted_2_bwd) {
        if (!exhausted_1_fwd) {
            wfa_it_1_fwd.extend([&](Offset q, Offset t, Offset ext) {
                bool exhausted = false;
                if (use_heuristics) {
                    auto q_e = q + ext;
                    auto t_e = t + ext;
                    if ((q_e == query_1.size() && t_e >= std::min(query_1.size(), target_1.size()))
                            || (t_e == target_1.size() && q_e >= std::min(query_1.size(), target_1.size()))) {
                        exhausted_1_fwd = true;
                        exhausted_2_fwd = true;
                        exhausted_1_bwd = true;
                        exhausted_2_bwd = true;
                        exhausted = true;
                    }
                }
                update_breakpoints<false, false, true>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff);
            }, breakpoints.min_p);

            if (wfa_it_1_fwd.get_p() >= breakpoints.min_p)
                break;
        }

        if (!exhausted_2_fwd) {
            wfa_it_2_fwd.extend([&](Offset q, Offset t, Offset ext) {
                bool exhausted = false;
                if (use_heuristics) {
                    auto q_e = q + ext;
                    auto t_e = t + ext;
                    if ((q_e == query_2.size() && t_e >= std::min(query_2.size(), target_2.size()))
                            || (t_e == target_2.size() && q_e >= std::min(query_2.size(), target_2.size()))) {
                        exhausted_1_fwd = true;
                        exhausted_2_fwd = true;
                        exhausted_1_bwd = true;
                        exhausted_2_bwd = true;
                        exhausted = true;
                    }
                }
                update_breakpoints<false, true, true>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff);
            }, breakpoints.min_p);

            if (wfa_it_2_fwd.get_p() >= breakpoints.min_p)
                break;
        }

        if (!exhausted_1_bwd) {
            wfa_it_1_bwd.extend([&](Offset q, Offset t, Offset ext) {
                bool exhausted = false;
                if (use_heuristics) {
                    auto q_e = q + ext;
                    auto t_e = t + ext;
                    if ((q_e == query_rc_1.size() && t_e >= std::min(target_1.size(), query_rc_1.size()))
                            || (t_e == target_1.size() && q_e >= std::min(target_1.size(), query_rc_1.size()))) {
                        exhausted_1_fwd = true;
                        exhausted_2_fwd = true;
                        exhausted_1_bwd = true;
                        exhausted_2_bwd = true;
                        exhausted = true;
                    }
                }
                update_breakpoints<true, true, false>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff);
            }, breakpoints.min_p);

            if (wfa_it_1_bwd.get_p() >= breakpoints.min_p)
                break;
        }

        if (!exhausted_2_bwd) {
            wfa_it_2_bwd.extend([&](Offset q, Offset t, Offset ext) {
                bool exhausted = false;
                if (use_heuristics) {
                    auto q_e = q + ext;
                    auto t_e = t + ext;
                    if ((q_e == query_rc_2.size() && t_e >= std::min(target_2.size(), query_rc_2.size()))
                            || (t_e == target_2.size() && q_e >= std::min(query_rc_2.size(), target_2.size()))) {
                        exhausted_1_fwd = true;
                        exhausted_2_fwd = true;
                        exhausted_1_bwd = true;
                        exhausted_2_bwd = true;
                        exhausted = true;
                    }
                }
                update_breakpoints<false, false, false>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff);
            }, breakpoints.min_p);

            if (wfa_it_2_bwd.get_p() >= breakpoints.min_p)
                break;
        }

        if (!exhausted_1_fwd)
            wfa_it_1_fwd.next();

        if (!exhausted_2_fwd)
            wfa_it_2_fwd.next();

        if (!exhausted_1_bwd)
            wfa_it_1_bwd.next();

        if (!exhausted_2_bwd)
            wfa_it_2_bwd.next();

        exhausted_1_fwd |= wfa_it_1_fwd.empty();
        exhausted_1_bwd |= wfa_it_1_bwd.empty();
        exhausted_2_fwd |= wfa_it_2_fwd.empty();
        exhausted_2_bwd |= wfa_it_2_bwd.empty();
    }

    return get_alignment_cigars(aligner,
                                score_model,
                                query_1,
                                query_rc_1,
                                target_1,
                                query_2,
                                query_rc_2,
                                target_2,
                                breakpoints,
                                heuristics_length_cutoff);
}