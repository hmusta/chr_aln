#include "chaining.hpp"

#include "repeat_aligner.hpp"

void print_chain(std::ostream &out,
                 std::string_view theader,
                 std::string_view target,
                 std::string_view qheader,
                 std::string_view query,
                 const std::vector<Ranges> &chain,
                 Score chain_score) {
    Diag target_diag = static_cast<Diag>(target.size()) - static_cast<Diag>(query.size());
    size_t begin_i = 0;
    auto front_it = chain.begin();
    if (front_it->size() == 0) {
        ++front_it;
        ++begin_i;
        assert(front_it != chain.end());
    }
    size_t end_i = chain.size();
    auto back_it = chain.rbegin();
    if (back_it->size() == 0) {
        ++back_it;
        --end_i;
        assert(back_it != chain.rend());
    }

    out << "Chain score:\t" << chain_score << " with "
        << (back_it.base() - front_it) << " anchors" << "\n"
        << "t: " << front_it->rbegin + 1 << "-" << back_it->rend
        << " / " << target.size() << "\t"
        << "q: " << front_it->qbegin + 1 << "-" << back_it->qend
        << " / " << query.size() << std::endl;


    for (size_t i = begin_i; i < end_i; ++i) {
        out << "i: " << i << "\t"
            << orientation_to_char(chain[i].rorientation) << "\t" << theader << ":"
            << chain[i].rbegin + 1 << "-" << chain[i].rend << "\t"
            << orientation_to_char(chain[i].qorientation) << "\t" << qheader << ":"
            << chain[i].qbegin + 1 << "-" << chain[i].qend << "\t"
            << chain[i].size() << "\t"
            << chain[i].diag() << " / " << target_diag << "\t"
            << "(" << chain[i].left_trim << "," << chain[i].right_trim << ")"
            << std::endl;
    }
}

std::pair<std::vector<Ranges>, Score> chain_ranges(std::string_view target,
                                                   std::string_view query,
                                                   std::string_view query_rc,
                                                   std::vector<Ranges>& ranges,
                                                   const ScoreModel &score_model,
                                                   const ChainScoreModel &chain_score_model,
                                                   bool chain_inversions,
                                                   size_t nthreads) {
    if (ranges.empty())
        return std::make_pair(ranges, 0);

    #ifndef NDEBUG
    for (const auto &range : ranges) {
        range.assert_valid(false);
    }
    #endif

    // add endpoints
    ranges.emplace_back(0, 0, false, 0, 0, false);
    ranges.emplace_back(target.size(), target.size(), false, query.size(), query.size(), false);

    auto get_overlap_score = [&chain_score_model](SOffset overlap, SOffset /* right_trim */) -> Score {
        return chain_score_model.get_match_score(overlap);
    };

    auto get_gap_score = [&chain_score_model](SOffset rend_a, SOffset rend_b, SOffset qend_a,
                                              SOffset qend_b, SOffset overlap = 0, bool is_end = false) -> Score {
        assert(rend_a <= rend_b);
        assert(qend_a <= qend_b);
        SOffset dist = rend_b - rend_a;
        SOffset qdist = qend_b - qend_a;

        Offset gap_edits = std::abs(dist - qdist);
        return chain_score_model.get_gap_score(gap_edits);
    };

    auto get_mismatch_score
            = [&chain_score_model](SOffset rend_a, SOffset rend_b, SOffset qend_a,
                                SOffset qend_b, SOffset overlap = 0, bool is_end = false) -> Score {
        assert(rend_a <= rend_b);
        assert(qend_a <= qend_b);
        SOffset dist = rend_b - rend_a;
        SOffset qdist = qend_b - qend_a;

        SOffset min_dist = std::min(dist, qdist);
        assert(min_dist >= overlap);

        return chain_score_model.get_mismatch_score(min_dist);
    };

    auto get_inv_score = [&chain_score_model]<bool open>(SOffset qend_a, SOffset qend_b) -> Score {
        assert(qend_a <= qend_b);
        SOffset inv_len = qend_b - qend_a;
        Score inv_score = inv_len * chain_score_model.inv_ext_s;

        if constexpr (open) {
            inv_score += chain_score_model.inv_open_s;
        }

        assert(inv_score <= 0);
        return inv_score;
    };

    std::vector<ChainTableElem> dp_table;
    dp_table.reserve(ranges.size());

    for (auto it = ranges.begin(); it != ranges.end(); ++it) {
        Score chain_score = (it->rend == 0) ? 0 : ScoreModel::ninf_s;
        dp_table.emplace_back(chain_score,
                              !it->qorientation ? ScoreModel::ninf_s : chain_score, // inverse chain score
                              it);
    }
    assert(dp_table.size() == ranges.size());

    std::sort(dp_table.begin(), dp_table.end());

    auto update = [&](size_t i, size_t j) -> uint8_t {
        dp_table[i].assert_valid();
        const auto& [chain_score_a, chain_score_inv_a, it_a, last_a, last_inv_a,
                     start_inv_a, last_op_a, last_op_inv_a] = dp_table[i];

        it_a->assert_valid();
        const auto& [rbegin_a, rend_a, rorientation_a, qbegin_a, qend_a, qorientation_a,
                     left_trim_a, right_trim_a] = *it_a;

        SOffset mum_length_a = it_a->size();

        dp_table[j].assert_valid();
        auto& [chain_score_b, chain_score_inv_b, it_b, last_b, last_inv_b, start_inv_b,
               last_op_b, last_op_inv_b] = dp_table[j];

        it_b->assert_valid();
        const auto& [rbegin_b, rend_b, rorientation_b, qbegin_b, qend_b, qorientation_b,
                     left_trim_b, right_trim_b] = *it_b;

        if (!chain_inversions && (qorientation_a || qorientation_b))
            return 0;

        SOffset mum_length_b = it_b->size();

        assert(rend_a <= rend_b);
        if (rend_a == rend_b && rbegin_a == rbegin_b && qbegin_a == qbegin_b && qend_a == qend_b)
            return 0;

        if (mum_length_a > 0 && mum_length_b > 0 && (rbegin_a >= rbegin_b || rend_a >= rend_b))
            return 0;

        if (!qorientation_a) { // forward to forward or inversion
            // make sure query is advanced
            if (mum_length_a > 0 && mum_length_b > 0 && (qend_a >= qend_b || qbegin_a >= qbegin_b))
                return 0;

            SOffset o1 = rend_b - std::max(rbegin_b, rend_a);
            SOffset o2 = qend_b - std::max(qbegin_b, qend_a);

            SOffset overlap = std::min(o1, o2);
            assert(overlap >= 0);
            assert(mum_length_b == 0 || overlap > 0);

            Score gap_score = get_gap_score(rend_a, rend_b, qend_a, qend_b, overlap,
                                              mum_length_a == 0 || mum_length_b == 0)
                                    + get_mismatch_score(rend_a, rend_b, qend_a, qend_b, overlap,
                                              mum_length_a == 0 || mum_length_b == 0);

            Score overlap_score_adj = get_overlap_score(overlap, right_trim_b);

            Score base_chain_score = chain_score_a;
            if (qorientation_b) {
                assert(chain_inversions);
                base_chain_score += get_inv_score.operator()<true>(qbegin_b, qend_b);
            }

            Score chain_score = base_chain_score + overlap_score_adj + gap_score;
            if (chain_score > chain_score_b) {
                chain_score_b = chain_score;
                last_b = i;

                if (qorientation_b) {
                    last_op_b = LastOp::INVERTED;
                    chain_score_inv_b = base_chain_score;
                    last_inv_b = i;
                    start_inv_b = j;
                    last_op_inv_b = LastOp::MATCH;
                } else {
                    last_op_b = LastOp::MATCH;
                }

                return 2;
            } else if (qorientation_b && base_chain_score > chain_score_inv_b) {
                last_op_b = LastOp::INVERTED;
                chain_score_inv_b = base_chain_score;
                last_inv_b = i;
                start_inv_b = j;
                last_op_inv_b = LastOp::MATCH;

                return 2;
            }

            return 1;
        }

        assert(chain_inversions);
        assert(mum_length_a > 0);
        if (start_inv_a == ChainTableElem::npos)
            return 0;

        assert(start_inv_a < dp_table.size());

        if (qorientation_b) { // inversion to the next inversion
            assert(mum_length_b > 0);

            // must advance query, so we need qbegin_b < qbegin_a and qend_b < qend_a
            if (qbegin_b >= qbegin_a || qend_b >= qend_a)
                return 0;

            assert(dp_table[start_inv_a].last_inv < dp_table.size());
            auto prev_it = dp_table[dp_table[start_inv_a].last_inv].it;
            assert(!prev_it->qorientation);
            SOffset prev_qend = prev_it->qend;
            SOffset prev_rend = prev_it->rend;

            size_t i_prev = prev_it - ranges.begin();
            size_t i_start = dp_table[start_inv_a].it - ranges.begin();
            size_t i_end = dp_table[j].it - ranges.begin();

            assert(!ranges[i_prev].qorientation);
            assert(ranges[i_start].qorientation);
            assert(ranges[i_end].qorientation);

            SOffset rbegin_start = ranges[i_start].rbegin;
            assert(prev_it->rbegin <= rbegin_start);
            assert(prev_it->size() == 0 || prev_it->rbegin < rbegin_start);
            assert(prev_rend < ranges[i_start].rend);
            assert(prev_it->qbegin <= ranges[i_start].qbegin);
            assert(prev_it->size() == 0 || prev_it->qbegin < ranges[i_start].qbegin);
            assert(prev_qend < ranges[i_start].qend);

            // if the last forward anchor touched or overlapped with opening RC anchor, there's nowhere to go
            if (prev_rend >= rbegin_start || prev_qend >= qbegin_a)
                return 0;

            auto [ql_1, query_w_1, query_rc_w_1, target_w_1, mum_length_1, q_begin, t_begin, eml_1]
                = extract_gap_seqs_switch<false>(target, query, query_rc,
                                                 ranges, i_prev, i_start, i_end);
            assert(query_rc_w_1.empty());

            // for reference
            // best_chain[i].trim_prefix(std::max<SOffset>(ti_1 + target_w_1.size() - best_chain[i].rbegin, 0), false);
            // best_chain[i].trim_suffix(std::max<SOffset>(best_chain[i].qend - qi_2, 0), false);
            // best_chain[j - 1].trim_prefix(std::max<SOffset>(qi_1 + query_w_1.size() - best_chain[j - 1].qbegin, 0), false);
            // best_chain[j - 1].trim_suffix(std::max<SOffset>(best_chain[j - 1].rend - ti_2, 0), false);

            Offset start_trim = std::max<SOffset>(t_begin + target_w_1.size() - ranges[i_start].rbegin, 0);
            Offset end_trim = std::max<SOffset>(q_begin + query_w_1.size() - ranges[i_end].qbegin, 0);

            if (i_start != i_end) {
                if (start_trim >= ranges[i_start].size())
                    return 0;

                if (end_trim >= ranges[i_end].size())
                    return 0;
            } else {
                if (std::max(start_trim, end_trim) >= ranges[i_start].size())
                    return 0;
            }

            SOffset o1 = rend_b - std::max(rbegin_b, rend_a);
            assert(o1 >= 0);
            assert(mum_length_b == 0 || o1 > 0);

            SOffset o2 = std::min(qbegin_a, qend_b) - qbegin_b;
            assert(o2 >= 0);
            assert(mum_length_b == 0 || o2 > 0);

            SOffset overlap = std::min(o1, o2);
            assert(overlap >= 0);
            assert(mum_length_b == 0 || overlap > 0);

            Score gap_score = get_gap_score(rend_a, rend_b, qbegin_b, qbegin_a, overlap,
                                                mum_length_a == 0 || mum_length_b == 0)
                                    + get_mismatch_score(rend_a, rend_b, qbegin_b, qbegin_a, overlap,
                                                mum_length_a == 0 || mum_length_b == 0);

            Score overlap_score_adj = get_overlap_score(overlap, right_trim_b);

            Score chain_score = chain_score_inv_a
                    + gap_score
                    + get_inv_score.operator()<false>(qbegin_b, qbegin_a)
                    + overlap_score_adj;

            SOffset t_end = t_begin + target_w_1.size() + eml_1;
            SOffset q_end = q_begin + query_w_1.size() + eml_1;

            Score close_gap_score =
                get_gap_score(t_begin, t_end, q_begin, q_end)
                    + get_mismatch_score(t_begin, t_end, q_begin, q_end);

            Score close_overlap_score_adj = get_overlap_score(eml_1, 0);

            if (chain_score + close_overlap_score_adj + close_gap_score > chain_score_b) {
                chain_score_b = chain_score + close_overlap_score_adj + close_gap_score;
                last_b = i;
                last_op_b = LastOp::INVERTED;

                chain_score_inv_b = chain_score;
                last_inv_b = i;
                start_inv_b = start_inv_a;
                last_op_inv_b = LastOp::INVERTED;

                return 2;
            } else if (chain_score > chain_score_inv_b) {
                chain_score_inv_b = chain_score;
                last_inv_b = i;
                start_inv_b = start_inv_a;
                last_op_inv_b = LastOp::INVERTED;
                return 2;
            }

            return 1;
        } else { // inversion to forward
            assert(dp_table[start_inv_a].last_inv < dp_table.size());

            auto prev_it = dp_table[dp_table[start_inv_a].last_inv].it;
            assert(!prev_it->qorientation);

            // SOffset prev_rend = prev_it->rend;
            SOffset prev_qend = prev_it->qend;
            if (prev_qend >= qbegin_b)
                return 0;

            auto start_it = dp_table[start_inv_a].it;
            const auto& [rbegin_start, rend_start, rrc_start, qbegin_start, qend_start, qrc_start, left_trim_start,
                         right_trim_start] = *start_it;
            assert(qbegin_a <= qbegin_start);
            assert(start_inv_a == i || qbegin_a < qbegin_start);
            assert(qend_a <= qend_start);
            assert(start_inv_a == i || qend_a < qend_start);
            assert(start_it->size() > 0);

            // make sure query is advanced
            if (mum_length_b > 0 && (qend_start >= qend_b || qbegin_start >= qbegin_b))
                return 0;

            size_t i_prev = prev_it - ranges.begin();
            size_t i_start = start_it - ranges.begin();
            size_t i_end = dp_table[i].it - ranges.begin();
            size_t i_next = dp_table[j].it - ranges.begin();
            assert(!ranges[i_prev].qorientation);
            assert(ranges[i_start].qorientation);
            assert(ranges[i_end].qorientation);
            assert(!ranges[i_next].qorientation);

            auto [ql_1, query_w_1, query_rc_w_1, target_w_1, mum_length_1, q_begin_open, t_begin_open, overlap]
                = extract_gap_seqs_switch<false>(target, query, query_rc,
                                                 ranges,
                                                 i_prev,
                                                 i_start,
                                                 i_end,
                                                 i_next);

            auto [ql_2, query_w_2, query_rc_w_2, target_w_2, mum_length_2, q_begin_close, t_begin_close, close_overlap]
                = extract_gap_seqs_switch<true>(target, query, query_rc,
                                                ranges,
                                                i_prev,
                                                i_start,
                                                i_end,
                                                i_next);

            // for reference
            // Offset start_prefix_trim = std::max<SOffset>(ti_1 + target_w_1.size() - best_chain[i].rbegin, 0);
            // Offset start_suffix_trim = std::max<SOffset>(best_chain[i].qend - qi_2, 0);
            // best_chain[i].trim_prefix(start_prefix_trim, false);
            // best_chain[i].trim_suffix(start_suffix_trim, false);

            // Offset end_prefix_trim = std::max<SOffset>(qi_1 + query_w_1.size() - best_chain[j - 1].qbegin, 0);
            // Offset end_suffix_trim = std::max<SOffset>(best_chain[j - 1].rend - ti_2, 0);
            // best_chain[j - 1].trim_prefix(end_prefix_trim, false);
            // best_chain[j - 1].trim_suffix(end_suffix_trim, false);

            Offset start_prefix_trim = std::max<SOffset>(t_begin_open + target_w_1.size() - ranges[i_start].rbegin, 0);
            Offset start_suffix_trim = std::max<SOffset>(ranges[i_start].qend - q_begin_close, 0);
            Offset end_prefix_trim = std::max<SOffset>(q_begin_open + query_w_1.size() - ranges[i_end].qbegin, 0);
            Offset end_suffix_trim = std::max<SOffset>(ranges[i_end].rend - t_begin_close, 0);

            if (i_start != i_end) {
                if (start_prefix_trim + start_suffix_trim >= ranges[i_start].size())
                    return 0;

                if (end_prefix_trim + end_suffix_trim >= ranges[i_end].size())
                    return 0;
            } else {
                if (std::max(start_prefix_trim, end_prefix_trim) + std::max(start_suffix_trim, end_suffix_trim) >= ranges[i_start].size())
                    return 0;
            }

            assert(overlap >= 0);
            assert(mum_length_b == 0 || overlap > 0);
            assert(close_overlap >= 0);
            assert(mum_length_b == 0 || close_overlap > 0);

            Score overlap_score_adj = get_overlap_score(overlap, 0);
            Score close_overlap_score_adj = get_overlap_score(close_overlap, 0);

            SOffset t_end_open = t_begin_open + target_w_1.size() + overlap;
            SOffset q_end_open_1 = q_begin_open + query_w_1.size() + overlap;
            SOffset q_end_open_2 = q_begin_open + query_rc_w_1.size() + overlap;

            SOffset t_end_close = t_begin_close + target_w_2.size() + close_overlap;
            SOffset q_end_close_1 = q_begin_close + query_w_2.size() + close_overlap;
            SOffset q_end_close_2 = q_begin_close + query_rc_w_2.size() + close_overlap;

            Score close_gap_score = std::max(
                get_gap_score(t_begin_open, t_end_open, q_begin_open, q_end_open_1)
                    + get_mismatch_score(t_begin_open, t_end_open, q_begin_open, q_end_open_1)
                    + get_gap_score(t_begin_close, t_end_close, q_begin_close, q_end_close_2)
                    + get_mismatch_score(t_begin_close, t_end_close, q_begin_close, q_end_close_2),
                get_gap_score(t_begin_open, t_end_open, q_begin_open, q_end_open_2)
                    + get_mismatch_score(t_begin_open, t_end_open, q_begin_open, q_end_open_2)
                    + get_gap_score(t_begin_close, t_end_close, q_begin_close, q_end_close_1)
                    + get_mismatch_score(t_begin_close, t_end_close, q_begin_close, q_end_close_1)
            );

            Score chain_score = chain_score_inv_a
                                    + close_gap_score
                                    + overlap_score_adj + close_overlap_score_adj;

            if (chain_score > chain_score_b) {
                chain_score_b = chain_score;
                last_b = i;
                last_op_b = LastOp::MATCH;
                return 2;
            }

            return 1;
        }

        return 0;
    };

    ProgressBar progress_bar(dp_table.size() * (dp_table.size() - 1) / 2, "Chaining");
    for (size_t i = 0; i < dp_table.size(); ++i) {
        #pragma omp parallel for num_threads(nthreads)
        for (size_t j = i + 1; j < dp_table.size(); ++j) {
            update(i, j);
        }

        progress_bar += dp_table.size() - (i + 1);
    }

    // backtrack
    size_t i = dp_table.size() - 1;
    Score best_score = dp_table[i].chain_score;

    if (best_score == ScoreModel::ninf_s)
        return {};

    assert(dp_table[i].last_op == LastOp::MATCH);

    std::vector<Ranges> ret_val;
    while (i < dp_table.size()) {
        ret_val.emplace_back(*dp_table[i].it);
        assert(dp_table[i].last_op != LastOp::CLIPPED || dp_table[i].last == ChainTableElem::npos);

        if (dp_table[i].last_op == LastOp::INVERTED) {
            i = dp_table[i].last_inv;
            assert(i < dp_table.size());
            while (dp_table[i].last_op_inv == LastOp::INVERTED) {
                ret_val.emplace_back(*dp_table[i].it);
                i = dp_table[i].last_inv;
                assert(i < dp_table.size());
            }
        } else {
            i = dp_table[i].last;
        }
    }

    assert(ranges.size() > 2);
    ranges.erase(ranges.end() - 2, ranges.end());

    if (std::all_of(ret_val.begin(), ret_val.end(), [](const auto& a) { return a.size() == 0; }))
        return {};

    assert(ret_val.size() > 2);

    std::reverse(ret_val.begin(), ret_val.end());

    return std::make_pair(std::move(ret_val), best_score);
}

std::vector<Ranges> rc_ranges(const std::vector<Ranges>& mummer_ranges, SOffset query_size) {
    std::vector<Ranges> rc_mummer_ranges;
    rc_mummer_ranges.reserve(mummer_ranges.size());
    for (const auto& [rbegin, rend, rrc, qbegin, qend, qrc, left_trim, right_trim] :
         mummer_ranges) {
        assert(query_size >= qbegin);
        assert(query_size >= qend);
        rc_mummer_ranges.emplace_back(rbegin, rend, rrc, query_size - qend,
                                      query_size - qbegin, !qrc, left_trim, right_trim);
    }

    return rc_mummer_ranges;
}

std::tuple<bool, std::string_view, std::string_view, SOffset, SOffset, SOffset, SOffset>
extract_gap_seqs_continue(std::string_view target,
                          std::string_view query,
                          std::string_view query_rc,
                          const std::vector<Ranges>& best_chain,
                          size_t i) {
    assert(i);
    SOffset qi = 0;
    SOffset qi_b = 0;
    SOffset ti = 0;
    SOffset mum_length = 0;
    bool qorientation_last = false;
    {
        const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation, left_trim,
                     right_trim] = best_chain[i - 1];
        qi_b = qbegin;
        qi = qend;
        ti = rend;
        qorientation_last = qorientation;
    }
    SOffset exact_match_length = max_offset;

    const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation, left_trim,
                 right_trim] = best_chain[i];
    assert(rend - rbegin == qend - qbegin);
    mum_length = rend - rbegin;

    assert(rbegin <= rend);
    assert(qbegin <= qend);
    assert(!rorientation);

    assert(ti <= rend);
    assert((qorientation && qorientation_last) || qi <= qend);

    SOffset tlen = rend - ti;
    std::string_view target_w(target.data() + ti, tlen);
    SOffset qlen = 0;

    assert(qorientation == qorientation_last);

    // TODO: add a case for switching orientation from inv to forward back to inv
    assert((qend >= qi) == !qorientation);

    // keep the same orientation
    assert(!qorientation || qi_b >= qbegin);
    qlen = !qorientation ? qend - qi : qi_b - qbegin;

    exact_match_length = std::min({ qlen, tlen, mum_length });
    qlen -= exact_match_length;
    tlen -= exact_match_length;

    // if qorientation, we want to go from qbegin to qi_b in query, so size - qi_b to size - qbegin in query_rc
    std::string_view query_w = !qorientation
            ? std::string_view(query.data() + qi, qlen)
            : std::string_view(query_rc.data() + query_rc.size() - qi_b, qlen);
    target_w.remove_suffix(exact_match_length);
    assert(static_cast<SOffset>(target_w.size()) == tlen);
    assert(static_cast<SOffset>(query_w.size()) == qlen);

    return std::make_tuple(qorientation_last, std::move(query_w), std::move(target_w),
                           mum_length, !qorientation ? qi : qi_b, ti, exact_match_length);
}

void reseed_large_gaps(std::string_view target,
                       std::string_view query,
                       std::string_view query_rc,
                       const ScoreModel &score_model,
                       std::vector<Ranges>& best_chain,
                       std::vector<size_t>& inv_starts,
                       std::vector<size_t>& inv_ends,
                       SOffset k,
                       const std::function<bool(SOffset, SOffset)> &check_if_heuristics,
                       SOffset max_gap,
                       double exp_mismatch_frac_between_mum_bp,
                       bool check_inversions,
                       size_t nthreads) {
    size_t old_chain_size = best_chain.size();
    for (size_t i = 1; i < best_chain.size(); ++i) {
        assert(inv_starts.size() == best_chain.size());
        assert(inv_ends.size() == best_chain.size());
        bool qorientation_last = best_chain[i - 1].qorientation;
        bool qorientation = best_chain[i].qorientation;

        // skip inversion close and inversion extensions since they're already handled by the inversion open case below
        if (qorientation_last)
            continue;

        if (!qorientation) { // staying in the same orientation
            auto [qorientation_last, query_w, target_w, mum_length, qi, ti, exact_match_length]
                    = extract_gap_seqs_continue(target, query, query_rc, best_chain, i);

            if (!check_if_heuristics(query_w.size(), target_w.size()))
                continue;

            SOffset tlen = target_w.size();
            SOffset qlen = query_w.size();

            if (tlen < k || qlen < k)
                continue;

            assert(qi + query_w.size() <= query_rc.size());
            std::string_view query_rc_w(query_rc.data() + query_rc.size() - qi - query_w.size(),
                                        query_w.size());
            assert(query_rc_w == reverse_complement(query_w));

            std::vector<Ranges> mums;
            size_t num_fw_seeds = 0;
            size_t num_rc_seeds = 0;
            call_mums(query_w, query_rc_w, target_w, k, [&](Ranges&& mum) {
                num_fw_seeds += !mum.qorientation;
                num_rc_seeds += mum.qorientation;
                mums.emplace_back(std::move(mum));
            });
            std::cout << "Continue i: " << i
                    << " reseeding found "
                    << num_fw_seeds << " fw seeds and "
                    << num_rc_seeds << " rc seeds in "
                    << ti+1 << "-" << ti + target_w.size() << "(" << target_w.size() << ")\t"
                    << qi+1 << "-" << qi + query_w.size() << "(" << query_w.size() << ")\n";

            if (mums.empty())
                continue;

            ChainScoreModel local_chain_score_model(
                k,
                score_model.gap_switch,
                exp_mismatch_frac_between_mum_bp,
                max_gap
            );

            auto [local_chain, local_chain_score] = chain_ranges(
                target_w, query_w, query_rc_w, mums, score_model, local_chain_score_model,
                check_inversions,
                std::min<size_t>(nthreads, mums.size() * (mums.size() - 1))
            );

            if (local_chain.empty()) {
                std::sort(mums.begin(), mums.end());
                std::cout << "Continue i: " << i << " no chain found.\n";
                std::cout << "Has large gap: " << has_large_gap(query_w, target_w) << "\n";
                continue;
            }

            assert(local_chain.size() > 2);

            auto jt_b = local_chain.begin() + 1;
            auto jt_e = local_chain.end() - 1;
            assert(jt_b < jt_e);

            size_t new_fw_seeds = 0;
            size_t new_rc_seeds = 0;
            std::for_each(jt_b, jt_e, [&](Ranges& mum) {
                new_fw_seeds += !mum.qorientation;
                new_rc_seeds += mum.qorientation;
                mum.shift_start(ti, qi);
                assert(mum.check_equal(target, query, query_rc));
            });

            assert(best_chain[i - 1].rend < jt_b->rend);
            assert((jt_e - 1)->rend < best_chain[i].rend);

            size_t new_anchor_count = std::distance(jt_b, jt_e);

            Offset new_rbegin = jt_b->rbegin;
            Offset new_rend = (jt_e - 1)->rend;
            assert(new_rbegin < new_rend);
            size_t new_ref_coverage = new_rend - new_rbegin;

            Offset new_qbegin = std::min_element(jt_b, jt_e, [](const auto &a, const auto &b) { return a.qbegin < b.qbegin; })->qbegin;
            Offset new_qend = std::max_element(jt_b, jt_e, [](const auto &a, const auto &b) { return a.qend < b.qend; })->qend;

            assert(new_qbegin < new_qend);
            size_t new_qry_coverage = new_qend - new_qbegin;

            std::cout << "Continue i: " << i << " "
                    << "(" << new_fw_seeds << " fwd and " << new_rc_seeds << " rc). "
                    << "rechaining locally: "
                    << new_anchor_count << " anchors. "
                    << "Ref coverage: "
                    << new_rbegin+1 << "-" << new_rend << " " << new_ref_coverage << " / " << target_w.size() << " "
                    << "Qry coverage: "
                    << new_qbegin+1 << "-" << new_qend << " " << new_qry_coverage << " / " << query_w.size()
                    << "\n";

            best_chain.insert(best_chain.begin() + i, jt_b, jt_e);

            std::tie(inv_starts, inv_ends) = compute_invs(target, query, query_rc, best_chain);

            --i;

            assert(best_chain[i + 1] == *jt_b);
        } else { // inversion open
            assert(!qorientation_last && qorientation);

            auto [ql_1, query_w_1, query_rc_w_1, target_w_1, mum_length_1, qi_1, ti_1, eml_1]
                = extract_gap_seqs_switch<false>(target, query, query_rc,
                                                best_chain,
                                                inv_starts[i] - 1,
                                                inv_starts[i],
                                                inv_ends[i],
                                                inv_ends[i] + 1);
            assert(ql_1 != qorientation);

            auto [ql_2, query_w_2, query_rc_w_2, target_w_2, mum_length_2, qi_2, ti_2, eml_2]
                = extract_gap_seqs_switch<true>(target, query, query_rc,
                                                best_chain,
                                                inv_starts[i] - 1,
                                                inv_starts[i],
                                                inv_ends[i],
                                                inv_ends[i] + 1);
            assert(ql_2 == qorientation);

            assert(query_w_1.size() == query_w_2.size());
            assert(query_rc_w_1.size() == query_rc_w_2.size());

            if (!check_if_heuristics(
                        std::max(query_w_1.size(), query_rc_w_1.size()) + std::max(query_w_2.size(), query_rc_w_2.size()),
                        target_w_1.size() + target_w_2.size()
                    )) {
                continue;
            }

            SOffset prev_mum_length = 0;
            SOffset next_mum_length = 0;
            ssize_t prev_offset = 0;
            ssize_t next_offset = 0;
            assert(ti_1 >= prev_mum_length);
            assert(qi_1 >= prev_mum_length);
            ti_1 -= prev_mum_length;
            qi_1 -= prev_mum_length;

            assert(target_w_1.data() <= target_w_2.data() + target_w_2.size());
            SOffset tlen = target_w_2.data() + target_w_2.size() - target_w_1.data() + prev_mum_length + next_mum_length;
            std::string_view target_w(target_w_1.data() - prev_mum_length, tlen);

            assert(query_w_1.data() <= query_rc_w_2.data() + query_rc_w_2.size());
            assert(query_rc_w_1.data() <= query_w_2.data() + query_w_2.size());

            SOffset qlen = query_rc_w_2.data() + query_rc_w_2.size() - query_w_1.data();
            assert(qlen == query_w_2.data() + query_w_2.size() - query_rc_w_1.data());
            qlen += prev_mum_length + next_mum_length;

            std::string_view query_w(query_w_1.data() - prev_mum_length, qlen);
            assert(std::string_view(target_w.data(), prev_mum_length)
                    == std::string_view(query_w.data(), prev_mum_length));
            assert(std::string_view(target_w.data() + tlen - next_mum_length, next_mum_length)
                    == std::string_view(query_w.data() + qlen - next_mum_length, next_mum_length));

            std::string_view query_rc_w(query_rc_w_1.data() - next_mum_length, qlen);
            assert(query_w == reverse_complement(query_rc_w));

            if (tlen < k || qlen < k)
                continue;

            // add mums
            std::vector<Ranges> mums;
            size_t num_fw_seeds = 0;
            size_t num_rc_seeds = 0;

            auto push_mum = [&](Ranges&& mum) {
                assert(mum.size());
                mums.emplace_back(std::move(mum));
            };

            call_mums(query_w, query_rc_w, target_w, k, push_mum);

            call_mums(query_w, query_rc_w, target_w, k, push_mum,
                      0, query_w_1.size(),
                      0, query_rc_w_1.size(),
                      0, target_w_1.size());

            SOffset prev_covered = best_chain[i].size();
            for (size_t j = i + 1; j <= inv_ends[i]; ++j) {
                assert(best_chain[j - 1].qorientation);
                assert(best_chain[j].qorientation);
                auto [qorientation_last_j, query_w_j, target_w_j, mum_length_j, qi_j, ti_j, exact_match_length_j]
                        = extract_gap_seqs_continue(target, query, query_rc, best_chain, j);
                assert(qorientation_last_j);
                assert(best_chain[j].rend > best_chain[j - 1].rend);
                prev_covered += best_chain[j].rend - std::max(best_chain[j].rbegin, best_chain[j - 1].rend);
                if (!check_if_heuristics(query_w_j.size(), target_w_j.size()))
                    continue;

                assert(query_rc_w.begin() <= query_w_j.begin());
                assert(query_w_j.end() <= query_rc_w.end());
                assert(target_w.begin() <= target_w_j.begin());
                assert(target_w_j.end() <= target_w.end());
                call_mums(query_w, query_rc_w, target_w, k, push_mum,
                          0, 0,
                          query_w_j.data() - query_rc_w.data(), query_w_j.data() - query_rc_w.data() + query_w_j.size(),
                          target_w_j.data() - target_w.data(), target_w_j.data() - target_w.data() + target_w_j.size());
            }

            assert(query_w.begin() <= query_rc_w_2.begin());
            assert(query_rc_w_2.end() <= query_w.end());
            assert(query_rc_w.begin() <= query_w_2.begin());
            assert(query_w_2.end() <= query_rc_w.end());
            assert(target_w.begin() <= target_w_2.begin());
            assert(target_w_2.end() <= target_w.end());
            call_mums(query_w, query_rc_w, target_w, k, push_mum,
                      query_w.size() - query_rc_w_2.size(), query_w.size(),
                      query_rc_w.size() - query_w_2.size(), query_rc_w.size(),
                      target_w.size() - target_w_2.size(),
                      target_w.size());

            std::sort(mums.begin(), mums.end());
            mums.erase(std::unique(mums.begin(), mums.end()), mums.end());

            for (const auto &mum : mums) {
                num_fw_seeds += !mum.qorientation;
                num_rc_seeds += mum.qorientation;
            }

            std::cout << "Switch i: " << i << " reseeding found "
                    << num_fw_seeds << " fw seeds and "
                    << num_rc_seeds << " rc seeds in "
                    << ti_1+1 << "-" << ti_1 + tlen << "(" << tlen << ")\t"
                    << qi_1+1 << "-" << qi_1 + qlen << "(" << qlen << ")\n";

            if (mums.empty())
                continue;

            ChainScoreModel local_chain_score_model(
                k,
                // local_min_mum_length,
                score_model.gap_switch,
                exp_mismatch_frac_between_mum_bp,
                max_gap
            );

            auto [local_chain, local_chain_score] = chain_ranges(
                target_w, query_w, query_rc_w, mums, score_model, local_chain_score_model,
                check_inversions,
                std::min<size_t>(nthreads, mums.size() * (mums.size() - 1))
            );

            if (local_chain.empty()) {
                std::sort(mums.begin(), mums.end());
                std::cout << "Switch i: " << i << " no chain found.\n";
                std::cout << "Has large gap: " << has_large_gap(query_w, target_w) << "\n";
                continue;
            }

            assert(local_chain.size() > 2);

            auto it_b = best_chain.begin() + inv_starts[i] + prev_offset;
            auto it_e = best_chain.begin() + inv_ends[i] + 1 + next_offset;
            assert(it_b < it_e);
            assert(std::all_of(it_b, it_e, [](const auto &a) { return a.qorientation; }));

            size_t prev_anchor_count = std::distance(it_b, it_e);

            auto jt_b = local_chain.begin() + 1;
            auto jt_e = local_chain.end() - 1;
            assert(jt_b < jt_e);

            size_t new_anchor_count = std::distance(jt_b, jt_e);
            size_t new_fw_seeds = 0;
            size_t new_rc_seeds = 0;
            std::for_each(jt_b, jt_e, [&](Ranges& mum) {
                assert(mum.size());
                new_fw_seeds += !mum.qorientation;
                new_rc_seeds += mum.qorientation;
                mum.shift_start(ti_1, qi_1);
                assert(mum.check_equal(target, query, query_rc));
            });

            SOffset new_covered = local_chain[1].size();
            for (size_t j = 2; j < local_chain.size() - 1; ++j) {
                assert(local_chain[j].rend > local_chain[j - 1].rend);
                new_covered += local_chain[j].rend - std::max(local_chain[j].rbegin, local_chain[j - 1].rend);
            }

            std::cout << "Switch i: [" << it_b - best_chain.begin() << ":" << it_e - best_chain.begin() << ") "
                    << "(" << new_fw_seeds << " fwd and " << new_rc_seeds << " rc). "
                    << "rechaining locally: "
                    << prev_anchor_count << " -> " << new_anchor_count << " anchors. "
                    << "Coverage "
                    << prev_covered << " -> " << new_covered
                    << " / " << target_w.size() << ","
                    << query_w.size()
                    << "\n";

            if (new_covered <= prev_covered) {
                std::cout << "Switch i: " << i << " skipping.\n";
                continue;
            }

            if (std::distance(it_b, it_e) >= std::distance(jt_b, jt_e)) {
                std::copy(jt_b, jt_e, it_b);
                best_chain.erase(it_b + std::distance(jt_b, jt_e), it_e);
            } else {
                std::copy(jt_b, jt_b + std::distance(it_b, it_e), it_b);
                best_chain.insert(it_e, jt_b + std::distance(it_b, it_e), jt_e);
            }

            std::tie(inv_starts, inv_ends) = compute_invs(target, query, query_rc, best_chain);

            --i;
        }
    }

    std::cout << "Updated chain: " << old_chain_size << " -> " << best_chain.size() << "\n";
}

std::pair<Diag, Diag> compute_diag(const std::string& target,
                                   const std::string& query,
                                   const std::string& query_rc,
                                   const std::vector<Ranges>& best_chain,
                                   const std::vector<size_t>& inv_starts,
                                   const std::vector<size_t>& inv_ends) {
    // checking diag
    SOffset len_sum_r = 0;
    SOffset len_sum_q = 0;
    Diag max_diag = 0;

    auto update_diag = [&](SOffset exact_match_length, SOffset rlen, SOffset qlen) {
        assert(exact_match_length != max_offset);

        len_sum_r += rlen + exact_match_length;
        assert(len_sum_r <= static_cast<SOffset>(target.size()));

        len_sum_q += qlen + exact_match_length;
        assert(len_sum_q <= static_cast<SOffset>(query.size()));

        Diag diag = len_sum_r - len_sum_q;
        if (abs(diag) > abs(max_diag))
            max_diag = diag;
    };

    for (size_t i = 1; i < best_chain.size(); ++i) {
        bool qorientation_last = best_chain[i - 1].qorientation;
        bool qorientation = best_chain[i].qorientation;

        if (qorientation == qorientation_last) {
            auto [qol, query_w, target_w, mum_length, qi, ti, exact_match_length]
                    = extract_gap_seqs_continue(target, query, query_rc, best_chain, i);
            assert(qol == qorientation_last);
            assert(ti == len_sum_r);
            assert(target.c_str() + len_sum_r == target_w.data());

            assert(qorientation || qi == len_sum_q);
            assert(qorientation || query.c_str() + len_sum_q == query_w.data());

            update_diag(exact_match_length, target_w.size(), query_w.size());
        } else {
            auto [qol, query_w, query_rc_w, target_w, mum_length,
                  qi, ti, exact_match_length]
                    = !qorientation_last
                        ? extract_gap_seqs_switch<false>(
                                target, query, query_rc, best_chain,
                                inv_starts[i] - 1,
                                inv_starts[i],
                                inv_ends[i],
                                inv_ends[i] + 1
                            )
                        : extract_gap_seqs_switch<true>(
                                target, query, query_rc, best_chain,
                                inv_starts[i - 1] - 1,
                                inv_starts[i - 1],
                                inv_ends[i - 1],
                                inv_ends[i - 1] + 1
                            );
            assert(qol == qorientation_last);
            std::string_view query_count = !qorientation_last ? query_w : query_rc_w;

            assert(ti == len_sum_r);
            assert(target.c_str() + len_sum_r == target_w.data());

            assert(qorientation_last || qi == len_sum_q);
            assert(qorientation_last || query.c_str() + len_sum_q == query_w.data());

            update_diag(
                exact_match_length,
                target_w.size(),
                query_count.size()
            );
        }
    }

    Diag target_abs_diag
            = static_cast<Diag>(target.size()) - static_cast<Diag>(query.size());
    if (target_abs_diag < 0) {
        target_abs_diag *= -1;
        max_diag *= -1;
    }

    return std::make_pair(max_diag, target_abs_diag);
}

std::pair<std::vector<size_t>, std::vector<size_t>>
compute_invs(std::string_view target,
             std::string_view query,
             std::string_view query_rc,
             std::vector<Ranges>& best_chain) {
    std::vector<size_t> inv_ends(best_chain.size(), best_chain.size());
    for (size_t j = best_chain.size(); j > 0; --j) {
        size_t i = j - 1;
        const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation,
                        left_trim, right_trim] = best_chain[i];

        if (qorientation) {
            if (j < best_chain.size()) {
                const auto& [rbegin_next, rend_next, rorientation_next, qbegin_next,
                                qend_next, qorientation_next, left_trim_next,
                                right_trim_next] = best_chain[j];
                inv_ends[i] = !qorientation_next ? i : inv_ends[j];
            } else {
                inv_ends[i] = i;
            }
            assert(inv_ends[i] >= i);
            assert(inv_ends[i] + 1 < best_chain.size());
        }
    }

    std::vector<size_t> inv_starts(best_chain.size(), best_chain.size());
    for (size_t i = 0; i < best_chain.size(); ++i) {
        const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation,
                        left_trim, right_trim] = best_chain[i];
        if (qorientation) {
            if (i) {
                const auto& [rbegin_last, rend_last, rorientation_last, qbegin_last,
                                qend_last, qorientation_last, left_trim_last,
                                right_trim_last] = best_chain[i - 1];
                inv_starts[i] = !qorientation_last ? i : inv_starts[i - 1];
            } else {
                inv_starts[i] = i;
            }
            assert(inv_starts[i] <= i);
            assert(i <= inv_ends[i]);
            assert(inv_starts[i] > 0);
            assert(inv_starts[i] < best_chain.size());
        }
    }

    // trim ends
    for (size_t i = 1; i < best_chain.size(); ++i) {
        if (!best_chain[i - 1].qorientation && best_chain[i].qorientation) {
            assert(i == inv_starts[i]);

            size_t j = inv_ends[i] + 1;
            assert(j < best_chain.size());
            assert(best_chain[j - 1].qorientation);
            assert(!best_chain[j].qorientation);

            auto [ql_1, query_w_1, query_rc_w_1, target_w_1, mum_length_1, qi_1, ti_1, eml_1]
                = extract_gap_seqs_switch<false>(target, query, query_rc,
                                                 best_chain,
                                                 i - 1, i, j - 1, j);
            assert(!ql_1);

            assert(ti_1 == best_chain[i - 1].rend);
            assert(qi_1 == best_chain[i - 1].qend);

            auto [ql_2, query_w_2, query_rc_w_2, target_w_2, mum_length_2, qi_2, ti_2, eml_2]
                = extract_gap_seqs_switch<true>(target, query, query_rc,
                                                best_chain,
                                                i - 1, i, j - 1, j);
            assert(ql_2);
            assert(ti_2 + static_cast<SOffset>(target_w_2.size()) == best_chain[j].rbegin);
            assert(qi_2 + static_cast<SOffset>(query_rc_w_2.size()) == best_chain[j].qbegin);

            Offset start_prefix_trim = std::max<SOffset>(ti_1 + target_w_1.size() - best_chain[i].rbegin, 0);
            Offset start_suffix_trim = std::max<SOffset>(best_chain[i].qend - qi_2, 0);
            Offset end_prefix_trim = std::max<SOffset>(qi_1 + query_w_1.size() - best_chain[j - 1].qbegin, 0);
            Offset end_suffix_trim = std::max<SOffset>(best_chain[j - 1].rend - ti_2, 0);

            if (i != j - 1) {
                best_chain[i].trim_prefix(start_prefix_trim, false);
                best_chain[i].trim_suffix(start_suffix_trim, false);

                best_chain[j - 1].trim_prefix(end_prefix_trim, false);
                best_chain[j - 1].trim_suffix(end_suffix_trim, false);
            } else {
                best_chain[i].trim_prefix(std::max(start_prefix_trim, end_prefix_trim), false);
                best_chain[i].trim_suffix(std::max(start_suffix_trim, end_suffix_trim), false);
            }
        }
    }

    return std::make_pair(std::move(inv_starts), std::move(inv_ends));
}
