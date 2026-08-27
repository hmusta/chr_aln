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

    bool all_on_right(Offset t1_max, Offset q1_max, Offset t2_max, Offset q2_max) const {
        return q1 + q1_left_gap == 0 && t1 + t1_left_gap == t1_max
                        && q2 == q2_max && (t2_max > 0 || q1_max > 0 || q2_max > 0);
    }

    std::tuple<Score, Offset, Cigar, Cigar>
    align_all_right(wfa::WFAligner& aligner,
                    const ScoreModel &score_model,
                    std::string_view target_2,
                    std::string_view query_2,
                    std::string_view query_rc_2,
                    bool get_cigar,
                    SOffset heuristics_length_cutoff,
                    bool penalty_to_score) const {
        // TODO
        std::ignore = get_cigar;

        if (!t2_left_gap && !t2_right_gap && !q1_right_gap && !q2_right_gap) {
            std::string query_2_cat(query_2);
            query_2_cat += query_rc_2;
            auto [p, cigar] = get_alignment(
                aligner, score_model, query_2_cat, target_2,
                penalty_to_score,
                heuristics_length_cutoff
            );

            Offset target_fw = cigar_get_target_pos(cigar, query_2.size(), target_2.size(), query_2_cat.size()).second;

            auto [cigar2f, cigar2b] = cigar_split(
                cigar, target_fw, query_2.size(),
                target_2.size(), query_2_cat.size()
            );

            return std::make_tuple(
                p + (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * query_2.size(),
                target_fw,
                std::move(cigar2f), std::move(cigar2b)
            );
        } else {
            auto [p2f, cigar2f] = get_alignment(
                aligner, score_model, get_q2_fw(query_2), get_t2_fw(target_2),
                penalty_to_score,
                heuristics_length_cutoff
            );
            if (t2_left_gap) {
                cigar2f.push(TARGET_CONSUME_OP, t2_left_gap);
                p2f += penalty_to_score
                    ? score_model.get_gap_score(t2_left_gap)
                    : score_model.get_gap_penalty(t2_left_gap);
            }
            if (q1_right_gap) {
                cigar2f.push(QUERY_CONSUME_OP, q1_right_gap);
                if (penalty_to_score) {
                    p2f += score_model.get_gap_score(q1_right_gap)
                        + score_model.inv_ext_s * q1_right_gap;
                } else {
                    p2f += score_model.get_gap_penalty(q1_right_gap)
                        + score_model.inv_ext_p * q1_right_gap;
                }
            }
            p2f += (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * get_q2_fw(query_2).size();

            auto [p2b, cigar2b] = get_alignment(
                aligner, score_model, get_q2_bw(query_rc_2), get_t2_bw(target_2),
                penalty_to_score,
                heuristics_length_cutoff
            );

            if (t2_right_gap) {
                Cigar cigar_tmp(TARGET_CONSUME_OP, t2_right_gap);
                cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar2b.begin()), std::make_move_iterator(cigar2b.end()));
                cigar2b = std::move(cigar_tmp);

                p2b += penalty_to_score
                    ? score_model.get_gap_score(t2_right_gap)
                    : score_model.get_gap_penalty(t2_right_gap);
            }
            if (q2_right_gap) {
                Cigar cigar_tmp(QUERY_CONSUME_OP, q2_right_gap);
                cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar2b.begin()), std::make_move_iterator(cigar2b.end()));
                cigar2b = std::move(cigar_tmp);

                if (penalty_to_score) {
                    p2b += score_model.get_gap_score(q2_right_gap);
                } else {
                    p2b += score_model.get_gap_penalty(q2_right_gap);
                }
            }

            return std::make_tuple(
                p2f + p2b,
                get_t2_fw(target_2).size() + t2_left_gap,
                std::move(cigar2f), std::move(cigar2b)
            );
        }
    }

    bool all_on_left(Offset t1_max, Offset q1_max, Offset t2_max, Offset q2_max) const {
        return q1 + q1_left_gap == q1_max
                            && q2 + q2_right_gap == 0
                            && t2 + t2_left_gap == t2_max && (t1_max > 0 || q1_max > 0 || q2_max > 0);
    }

    std::tuple<Score, Offset, Cigar, Cigar>
    align_all_left(wfa::WFAligner& aligner,
                   const ScoreModel &score_model,
                   std::string_view target_1,
                   std::string_view query_1,
                   std::string_view query_rc_1,
                   bool get_cigar,
                   SOffset heuristics_length_cutoff,
                   bool penalty_to_score) const {
        // TODO
        std::ignore = get_cigar;

        if (!t1_left_gap && !t1_right_gap && !q1_left_gap && !q2_left_gap) {
            std::string query_1_cat(query_1);
            query_1_cat += query_rc_1;

            auto [p, cigar] = get_alignment(
                aligner, score_model, query_1_cat, target_1,
                penalty_to_score,
                heuristics_length_cutoff
            );

            Offset target_fw = cigar_get_target_pos(cigar, query_1.size(), target_1.size(), query_1_cat.size()).second;

            auto [cigar1f, cigar1b] = cigar_split(
                cigar, target_fw, query_1.size(),
                target_1.size(), query_1_cat.size()
            );

            return std::make_tuple(
                p + (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * query_rc_1.size(),
                target_fw,
                std::move(cigar1f), std::move(cigar1b)
            );
        } else {
            auto [p1f, cigar1f] = get_alignment(
                aligner, score_model, get_q1_fw(query_1), get_t1_fw(target_1),
                penalty_to_score,
                heuristics_length_cutoff
            );
            if (t1_left_gap) {
                cigar1f.push(TARGET_CONSUME_OP, t1_left_gap);
                p1f += penalty_to_score
                    ? score_model.get_gap_score(t1_left_gap)
                    : score_model.get_gap_penalty(t1_left_gap);
            }

            if (q1_left_gap) {
                cigar1f.push(QUERY_CONSUME_OP, q1_left_gap);
                if (penalty_to_score) {
                    p1f += score_model.get_gap_score(q1_left_gap);
                } else {
                    p1f += score_model.get_gap_penalty(q1_left_gap);
                }
            }
            auto [p1b, cigar1b] = get_alignment(
                aligner, score_model, get_q1_bw(query_rc_1), get_t1_bw(target_1),
                penalty_to_score,
                heuristics_length_cutoff
            );
            if (t1_right_gap) {
                Cigar cigar_tmp(TARGET_CONSUME_OP, t1_right_gap);
                cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar1b.begin()), std::make_move_iterator(cigar1b.end()));
                cigar1b = std::move(cigar_tmp);

                p1b += penalty_to_score
                    ? score_model.get_gap_score(t1_right_gap)
                    : score_model.get_gap_penalty(t1_right_gap);
            }
            if (q2_left_gap) {
                Cigar cigar_tmp(QUERY_CONSUME_OP, q2_left_gap);
                cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar1b.begin()), std::make_move_iterator(cigar1b.end()));
                cigar1b = std::move(cigar_tmp);

                if (penalty_to_score) {
                    p1b += score_model.get_gap_score(q2_left_gap)
                        + score_model.inv_ext_s * q2_left_gap;
                } else {
                    p1b += score_model.get_gap_penalty(q2_left_gap)
                        + score_model.inv_ext_p * q2_left_gap;
                }
            }

            p1b += (penalty_to_score ? score_model.inv_ext_s : score_model.inv_ext_p) * get_q1_bw(query_rc_1).size();

            return std::make_tuple(
                p1f + p1b,
                get_t1_fw(target_1).size() + t1_left_gap,
                std::move(cigar1f), std::move(cigar1b)
            );
        }
    }
};

inline std::ostream& operator<<(std::ostream& out, const Breakpoints &bp) {
    out << bp.min_p << "\t"
        << "t1: " << bp.t1 << " " << bp.t1_left_gap << "," << bp.t1_right_gap << "\t"
        << "q1: " << bp.q1 << " " << bp.q1_left_gap << "," << bp.q1_right_gap << "\t"
        << "t2: " << bp.t2 << " " << bp.t2_left_gap << "," << bp.t2_right_gap << "\t"
        << "q2: " << bp.q2 << " " << bp.q2_left_gap << "," << bp.q2_right_gap;
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
                        SOffset heuristics_length_cutoff,
                        bool check_frp) {
    assert(wfa_it_1_fwd_in.max_q() == wfa_it_2_fwd_in.max_q());
    assert(wfa_it_1_bwd_in.max_q() == wfa_it_2_bwd_in.max_q());
    assert(wfa_it_1_fwd_in.max_t() == wfa_it_1_bwd_in.max_t());
    assert(wfa_it_2_fwd_in.max_t() == wfa_it_2_bwd_in.max_t());

    const SOffset q1_max = wfa_it_1_fwd_in.max_q();
    const SOffset q2_max = wfa_it_1_bwd_in.max_q();
    const SOffset t1_max = wfa_it_1_fwd_in.max_t();
    const SOffset t2_max = wfa_it_2_fwd_in.max_t();

    auto update = [&aligner,&breakpoints,&score_model,q1_max,q2_max,t1_max,t2_max,
                   q,t,ext,check_frp,
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

        #ifdef NDEBUG
        std::ignore = q1p;
        std::ignore = q2p;
        std::ignore = t1p;
        std::ignore = t2p;
        #endif

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

        if (bp_new.all_on_right(t1_max, q1_max, t2_max, q2_max)) {
            // t1 will align to nothing, t2 will align to the concatenation of query_2 and query_rc_2
            if (t2_max == 0) {
                // penalty already covers query_rc_2, now we need to cover query_2
                update_p += score_model.get_gap_penalty(q1_max)
                            + score_model.inv_ext_p * q1_max;
            } else {
                auto [newp, target_fw, dummy_cigar_fw, dummy_cigar_bw] = bp_new.align_all_right(
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
        } else if (bp_new.all_on_left(t1_max, q1_max, t2_max, q2_max)) {
            // t1 will align to the concatenation of query_1 and query_rc_1, t2 will align to nothing
            if (t1_max == 0) {
                // penalty already covers query_1, now we need to cover query_rc_1
                update_p += score_model.get_gap_penalty(q2_max)
                            + score_model.inv_ext_p * q2_max;
            } else {
                auto [newp, target_fw, dummy_cigar_fw, dummy_cigar_bw] = bp_new.align_all_left(
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

        if (update_p >= breakpoints.min_p)
            return;

        breakpoints = bp_new;
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

    if (check_frp) {
        auto [q1, t1, p1, ext_1f] = wfa_it_1_fwd_in.get_frp();
        auto [q2p, t1p, p1p, ext_1b] = wfa_it_1_bwd_in.get_frp();
        auto [q1p, t2, p2, ext_2f] = wfa_it_2_fwd_in.get_frp();
        auto [q2, t2p, p2p, ext_2b] = wfa_it_2_bwd_in.get_frp();

        SOffset q1_gap = q1_max - q1 - q1p;
        SOffset q2_gap = q2_max - q2 - q2p;
        SOffset t1_gap = t1_max - t1 - t1p;
        SOffset t2_gap = t2_max - t2 - t2p;

        // Penalty cur_p;
        // if constexpr (!inv_open && !inv_ext && is_fwd) {
        //     cur_p = wfa_it_1_fwd_in.get_p();
        // } else if constexpr (inv_open && inv_ext && !is_fwd) {
        //     cur_p = wfa_it_1_bwd_in.get_p();
        // } else if constexpr (!inv_open && inv_ext && is_fwd) {
        //     cur_p = wfa_it_2_fwd_in.get_p();
        // } else {
        //     cur_p = wfa_it_2_bwd_in.get_p();
        // }

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
                continue;

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
            if (min_update_p_3 >= breakpoints.min_p)
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

std::tuple<Score, Cigar, Cigar, SOffset, SOffset, SOffset, SOffset,
           Score, Cigar, Cigar, SOffset, SOffset, SOffset, SOffset>
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
                                                 Score inv_ext_s = 0,
                                                 Diag min_k = min_diag,
                                                 Diag max_k = max_diag) -> std::pair<Score, Cigar> {
        auto [score, cigar] = get_alignment(
            aligner, score_model,
            query, target,
            true,
            heuristics_length_cutoff,
            min_k, max_k
        );
        score += inv_ext_s * query.size();
        return std::make_pair(score, std::move(cigar));
    };

    Offset r_consumed_1 = 0;
    Offset r_consumed_1_bw = 0;
    Offset q_consumed_1 = 0;
    Offset q_consumed_1_bw = 0;
    Offset r_consumed_2 = 0;
    Offset r_consumed_2_bw = 0;
    Offset q_consumed_2 = 0;
    Offset q_consumed_2_bw = 0;

    Cigar cigar_1;
    Cigar cigar_1_bw;
    Cigar cigar_2;
    Cigar cigar_2_bw;

    Score score_1_fw = 0;
    Score score_1_bw = 0;
    Score score_2_fw = 0;
    Score score_2_bw = 0;

    if (breakpoints.all_on_right(target_1.size(), query_1.size(), target_2.size(), query_rc_1.size())) {
        // t1 will align to nothing, t2 will align to the concatenation of query_2 and query_rc_2
        if (target_1.size()) {
            cigar_1 = std::to_string(target_1.size()) + TARGET_CONSUME_OP;
            r_consumed_1 = target_1.size();
            score_1_fw = score_model.get_gap_score(target_1.size());
        }

        std::tie(score_2_fw, r_consumed_2, cigar_2, cigar_2_bw) = breakpoints.align_all_right(
            aligner,
            score_model,
            target_2,
            query_2,
            query_rc_2,
            true,
            heuristics_length_cutoff,
            true
        );
        assert(r_consumed_2 <= target_2.size());
        r_consumed_2_bw = target_2.size() - r_consumed_2;
        q_consumed_2 = query_2.size();
        q_consumed_2_bw = query_rc_2.size();

        assert(r_consumed_1 + r_consumed_1_bw == target_1.size());
        assert(q_consumed_1 + q_consumed_2 == query_1.size());
        assert(q_consumed_1_bw + q_consumed_2_bw == query_rc_1.size());
    } else if (breakpoints.all_on_left(target_1.size(), query_1.size(), target_2.size(), query_rc_1.size())) {
        if (target_2.size()) {
            cigar_2 = std::to_string(target_2.size()) + TARGET_CONSUME_OP;
            r_consumed_2 = target_2.size();
            score_2_bw = score_model.get_gap_score(target_2.size());
        }
        std::tie(score_1_fw, r_consumed_1, cigar_1, cigar_1_bw) = breakpoints.align_all_left(
            aligner,
            score_model,
            target_1,
            query_1,
            query_rc_1,
            true,
            heuristics_length_cutoff,
            true
        );
        assert(r_consumed_1 <= target_1.size());
        r_consumed_1_bw = target_1.size() - r_consumed_1;
        q_consumed_1 = query_1.size();
        q_consumed_1_bw = query_rc_1.size();

        assert(r_consumed_2 + r_consumed_2_bw == target_2.size());
        assert(q_consumed_1 + q_consumed_2 == query_1.size());
        assert(q_consumed_1_bw + q_consumed_2_bw == query_rc_1.size());
    } else {
        r_consumed_1 = breakpoints.get_t1_fw(target_1).size();
        q_consumed_1 = breakpoints.get_q1_fw(query_1).size();
        std::tie(score_1_fw, cigar_1) = align_segment(
            breakpoints.get_q1_fw(query_1),
            breakpoints.get_t1_fw(target_1)
            // wfa_it_1_fwd.get_min_diag(),
            // wfa_it_1_fwd.get_max_diag()
        );

        if (breakpoints.t1_left_gap > 0) {
            cigar_1.push(TARGET_CONSUME_OP, breakpoints.t1_left_gap);
            score_1_fw += score_model.get_gap_score(breakpoints.t1_left_gap);
            r_consumed_1 += breakpoints.t1_left_gap;
        }

        if (breakpoints.q1_left_gap > 0) {
            cigar_1.push(QUERY_CONSUME_OP, breakpoints.q1_left_gap);
            score_1_fw += score_model.get_gap_score(breakpoints.q1_left_gap);
            q_consumed_1 += breakpoints.q1_left_gap;
        }

        assert(cigar_1
            == cigar_fix_n(cigar_1, target_1.substr(0, r_consumed_1), query_1.substr(0, q_consumed_1)));

        r_consumed_1_bw = breakpoints.get_t1_bw(target_1).size();
        q_consumed_1_bw = breakpoints.get_q1_bw(query_rc_1).size();
        std::tie(score_1_bw, cigar_1_bw) = align_segment(
            breakpoints.get_q1_bw(query_rc_1),
            breakpoints.get_t1_bw(target_1),
            score_model.inv_ext_s
        );

        if (breakpoints.t1_right_gap > 0) {
            Cigar cigar_tmp(TARGET_CONSUME_OP, breakpoints.t1_right_gap);
            cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar_1_bw.begin()), std::make_move_iterator(cigar_1_bw.end()));
            cigar_1_bw = std::move(cigar_tmp);
            score_1_bw += score_model.get_gap_score(breakpoints.t1_right_gap);
            r_consumed_1_bw += breakpoints.t1_right_gap;
        }

        if (breakpoints.q2_left_gap > 0) {
            Cigar cigar_tmp(QUERY_CONSUME_OP, breakpoints.q2_left_gap);
            cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar_1_bw.begin()), std::make_move_iterator(cigar_1_bw.end()));
            cigar_1_bw = std::move(cigar_tmp);
            score_1_bw += score_model.get_gap_score(breakpoints.q2_left_gap) + score_model.inv_ext_s * breakpoints.q2_left_gap;
            q_consumed_1_bw += breakpoints.q2_left_gap;
        }

        assert(cigar_1_bw
            == cigar_fix_n(cigar_1_bw, target_1.substr(r_consumed_1), query_rc_1.substr(query_rc_1.size() - q_consumed_1_bw)));

        r_consumed_2 = breakpoints.get_t2_fw(target_2).size();
        q_consumed_2 = breakpoints.get_q2_fw(query_2).size();
        std::tie(score_2_fw, cigar_2) = align_segment(
            breakpoints.get_q2_fw(query_2),
            breakpoints.get_t2_fw(target_2),
            score_model.inv_ext_s
            // wfa_it_2_fwd.get_min_diag(),
            // wfa_it_2_fwd.get_max_diag()
        );

        if (breakpoints.t2_left_gap > 0) {
            cigar_2.push(TARGET_CONSUME_OP, breakpoints.t2_left_gap);
            score_2_fw += score_model.get_gap_score(breakpoints.t2_left_gap);
            r_consumed_2 += breakpoints.t2_left_gap;
        }

        if (breakpoints.q1_right_gap > 0) {
            cigar_2.push(QUERY_CONSUME_OP, breakpoints.q1_right_gap);
            score_2_fw += score_model.get_gap_score(breakpoints.q1_right_gap) + score_model.inv_ext_s * breakpoints.q1_right_gap;
            q_consumed_2 += breakpoints.q1_right_gap;
        }

        assert(cigar_2
            == cigar_fix_n(cigar_2, target_2.substr(0, r_consumed_2), query_2.substr(0, q_consumed_2)));

        r_consumed_2_bw = breakpoints.get_t2_bw(target_2).size();
        q_consumed_2_bw = breakpoints.get_q2_bw(query_rc_2).size();
        std::tie(score_2_bw, cigar_2_bw) = align_segment(
            breakpoints.get_q2_bw(query_rc_2),
            breakpoints.get_t2_bw(target_2)
        );

        if (breakpoints.t2_right_gap > 0) {
            Cigar cigar_tmp(TARGET_CONSUME_OP, breakpoints.t2_right_gap);
            cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar_2_bw.begin()), std::make_move_iterator(cigar_2_bw.end()));
            cigar_2_bw = std::move(cigar_tmp);

            score_2_bw += score_model.get_gap_score(breakpoints.t2_right_gap);
            r_consumed_2_bw += breakpoints.t2_right_gap;
        }

        if (breakpoints.q2_right_gap > 0) {
            Cigar cigar_tmp(QUERY_CONSUME_OP, breakpoints.q2_right_gap);
            cigar_tmp.insert(cigar_tmp.end(), std::make_move_iterator(cigar_2_bw.begin()), std::make_move_iterator(cigar_2_bw.end()));
            cigar_2_bw = std::move(cigar_tmp);

            score_2_bw += score_model.get_gap_score(breakpoints.q2_right_gap);
            q_consumed_2_bw += breakpoints.q2_right_gap;
        }

        assert(cigar_2_bw
            == cigar_fix_n(cigar_2_bw, target_2.substr(r_consumed_2), query_rc_2.substr(query_rc_2.size() - q_consumed_2_bw)));

        assert(r_consumed_1 + r_consumed_1_bw == target_1.size());
        assert(r_consumed_2 + r_consumed_2_bw == target_2.size());
        assert(q_consumed_1 + q_consumed_2 == query_1.size());
        assert(q_consumed_1_bw + q_consumed_2_bw == query_rc_1.size());
    }

    return std::make_tuple(
        score_1_fw + score_model.inv_open_s + score_1_bw,
        std::move(cigar_1), std::move(cigar_1_bw),
        r_consumed_1 + r_consumed_1_bw,
        q_consumed_1 + q_consumed_1_bw,
        q_consumed_1_bw,
        r_consumed_1_bw,
        score_2_fw + score_2_bw,
        std::move(cigar_2), std::move(cigar_2_bw),
        r_consumed_2 + r_consumed_2_bw,
        q_consumed_2 + q_consumed_2_bw,
        q_consumed_2,
        r_consumed_2
    );
}

std::tuple<Score, Cigar, SOffset, SOffset, SOffset, SOffset,
           Score, Cigar, SOffset, SOffset, SOffset, SOffset,
           SOffset, SOffset>
run_alignment(wfa::WFAligner& aligner,
              const ScoreModel &score_model,
              std::string_view query,
              std::string_view query_rc,
              std::string_view target,
              std::string_view query_1,
              std::string_view query_rc_1,
              std::string_view target_1,
              std::string_view query_2,
              std::string_view query_rc_2,
              std::string_view target_2,
              SOffset heuristics_length_cutoff,
              SOffset min_wavefront_length,
              SOffset max_distance_threshold) {
    assert(is_reverse_complement(query_1, query_2));
    assert(is_reverse_complement(query_rc_1, query_rc_2));

    assert(target_1.data() >= target.data());
    assert(target_2.data() >= target_1.data() + target_1.size());
    assert(target.data() + target.size() >= target_2.data() + target_2.size());

    assert(query_1.data() >= query.data());
    assert(query_rc_2.data() >= query_1.data() + query_1.size());
    assert(query.data() + query.size() >= query_rc_2.data() + query_rc_2.size());

    assert(query_rc_1.data() >= query_rc.data());
    assert(query_2.data() >= query_rc_1.data() + query_rc_1.size());
    assert(query_rc.data() + query_rc.size() >= query_2.data() + query_2.size());
    #ifdef NDEBUG
    std::ignore = query_rc;
    #endif

    Offset qi_1 = query_1.data() - query.data();
    Offset ti_1 = target_1.data() - target.data();
    Offset ti_2 = target_2.data() - target.data();

    Breakpoints breakpoints {
        .min_p = ScoreModel::inf_p,
        .t1 = target_1.size() + 1,
        .t2 = target_2.size() + 1,
        .q1 = query_1.size() + 1,
        .q2 = query_rc_1.size() + 1
    };

    const bool use_heuristics = (min_wavefront_length < max_diag_width)
                                || (max_distance_threshold < max_offset);

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

    update_breakpoints<false, false, true>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, 0, 0, 0, false, heuristics_length_cutoff, true);

    bool exhausted_1_fwd = target_1.empty() || query_1.empty();
    bool exhausted_1_bwd = target_1.empty() || query_rc_1.empty();
    bool exhausted_2_fwd = target_2.empty() || query_2.empty();
    bool exhausted_2_bwd = target_2.empty() || query_rc_2.empty();

    UDiag max_antidiag_1_fwd = wfa_it_1_fwd.get_max_antidiag();
    UDiag max_antidiag_1_bwd = wfa_it_1_bwd.get_max_antidiag();
    UDiag max_antidiag_2_fwd = wfa_it_2_fwd.get_max_antidiag();
    UDiag max_antidiag_2_bwd = wfa_it_2_bwd.get_max_antidiag();

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

                bool check_frp = false;
                if (wfa_it_1_fwd.get_max_antidiag() > max_antidiag_1_fwd) {
                    max_antidiag_1_fwd = wfa_it_1_fwd.get_max_antidiag();
                    check_frp = true;
                }

                update_breakpoints<false, false, true>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff, check_frp);
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

                bool check_frp = false;
                if (wfa_it_2_fwd.get_max_antidiag() > max_antidiag_2_fwd) {
                    max_antidiag_2_fwd = wfa_it_2_fwd.get_max_antidiag();
                    check_frp = true;
                }

                update_breakpoints<false, true, true>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff, check_frp);
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

                bool check_frp = false;
                if (wfa_it_1_bwd.get_max_antidiag() > max_antidiag_1_bwd) {
                    max_antidiag_1_bwd = wfa_it_1_bwd.get_max_antidiag();
                    check_frp = true;
                }

                update_breakpoints<true, true, false>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff, check_frp);
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

                bool check_frp = false;
                if (wfa_it_2_bwd.get_max_antidiag() > max_antidiag_2_bwd) {
                    max_antidiag_2_bwd = wfa_it_2_bwd.get_max_antidiag();
                    check_frp = true;
                }

                update_breakpoints<false, false, false>(aligner, score_model, breakpoints, wfa_it_1_fwd, wfa_it_1_bwd, wfa_it_2_fwd, wfa_it_2_bwd, q, t, ext, exhausted, heuristics_length_cutoff, check_frp);
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

    auto [score_1, cigar_1, cigar_1_bw,
          r_consumed_1, q_consumed_1,
          inverted_1, inverted_r_1,
          score_2, cigar_2, cigar_2_bw,
          r_consumed_2, q_consumed_2,
          inverted_2, inverted_r_2] = get_alignment_cigars(aligner,
                                                            score_model,
                                                            query_1,
                                                            query_rc_1,
                                                            target_1,
                                                            query_2,
                                                            query_rc_2,
                                                            target_2,
                                                            breakpoints,
                                                            heuristics_length_cutoff);

    Offset q_inv_begin = qi_1 + q_consumed_1 - inverted_1;
    Offset q_inv_end = query_rc_2.data() - query.data() + query_rc_2.size() - (q_consumed_2 - inverted_2);

    assert(q_inv_begin < q_inv_end);
    SOffset inv_len = q_inv_end - q_inv_begin;

    Offset r_inv_begin = ti_1 + r_consumed_1 - inverted_r_1;
    Offset r_inv_end = ti_2 + target_2.size() - (r_consumed_2 - inverted_r_2);

    assert(r_inv_begin < r_inv_end);
    SOffset inv_len_r = r_inv_end - r_inv_begin;

    const char *q_inv_data = query_rc_1.data() + query_rc_1.size() - inverted_1;

    assert(target_1.data() + r_consumed_1 - inverted_r_1 == target.data() + r_inv_begin);
    assert(query_1.data() + q_consumed_1 - inverted_1 == query.data() + q_inv_begin);

    assert(target_2.data() + target_2.size() - (r_consumed_2 - inverted_r_2) == target.data() + r_inv_end);
    assert(query_rc_2.data() + query_rc_2.size() - (q_consumed_2 - inverted_2) == query.data() + q_inv_end);

    std::string_view q_inv_window(q_inv_data, inv_len);
    assert(q_inv_data + inv_len == query_2.data() + inverted_2);
    assert(is_reverse_complement(q_inv_window, std::string_view(query.data() + q_inv_begin, inv_len)));

    cigar_1.push(INV_OP, inv_len);
    cigar_1.insert(cigar_1.end(), std::make_move_iterator(cigar_1_bw.begin()), std::make_move_iterator(cigar_1_bw.end()));
    cigar_2.insert(cigar_2.end(), std::make_move_iterator(cigar_2_bw.begin()), std::make_move_iterator(cigar_2_bw.end()));

    return std::make_tuple(
        score_1, std::move(cigar_1), r_consumed_1, q_consumed_1, inverted_1, inverted_r_1,
        score_2, std::move(cigar_2), r_consumed_2, q_consumed_2, inverted_2, inverted_r_2,
        inv_len, inv_len_r
    );
}