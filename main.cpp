#include "chaining.hpp"
#include "cigar.hpp"
#include "helpers.hpp"
#include "readers.hpp"
#include "repeat_aligner.hpp"
#include "wfa_main.hpp"
#include "wfa_switch.hpp"

#ifdef PROFILING
#include "profiler.hpp"
#endif

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <cassert>
#include <cstdint>
#include <omp.h>

#include <progress_bar.hpp>

int main(int argc, char** argv) {
    if (argc <= 6) {
        std::cerr << "Usage: chr_aln REF.fa QRY.fa MUMMER.out NTHREADS CHECK_INVERSIONS CHAIN_OUT.out [CHAIN_QRY_RC] [ACCURATE_MODE]"
                  << std::endl;
        return 1;
    }

    #ifdef PROFILING
    auto profile_flusher = start_profiler_flusher();
    #endif

    std::cout << "Parsing sequences\n";
    assert(argc > 1);
    std::ifstream fint(argv[1]);
    auto [target, theader] = read_fasta(fint);
    fint.close();

    assert(argc > 2);
    std::ifstream finq(argv[2]);
    auto [query, qheader] = read_fasta(finq);
    finq.close();

    Diag target_diag = static_cast<Diag>(target.size()) - static_cast<Diag>(query.size());

    ScoreModel score_model(1, -9, -16, -2, -41, -1, -41, 0);

    double exp_mismatch_frac_between_mum_bp = 0.01;
    SOffset max_gap = max_offset;

    SOffset k = 31;

    wfa::WFAligner::MemoryModel heuristics_model = wfa::WFAligner::MemoryLow;
    wfa::WFAligner::MemoryModel main_model = wfa::WFAligner::MemoryUltralow;

    std::cout << score_model << std::endl;

    // TODO: handle later
    assert(score_model.inv_ext_p == 0);
    assert(score_model.inv_ext_s == 0);
    if (score_model.inv_ext_s != 0 || score_model.inv_ext_p != 0) {
        throw std::runtime_error("Non-zero inversion extension not supported yet");
    }

    /* Setting heuristics */
    // SOffset max_dist = std::numeric_limits<SOffset>::max();
    SOffset max_dist = 50;
    SOffset min_wavefront_length = 10;
    UDiag max_diag_diff = std::numeric_limits<UDiag>::max();

    size_t long_seq_length_cutoff = 100000;
    // size_t long_seq_length_cutoff = std::numeric_limits<SOffset>::max();

    SOffset heuristics_length_cutoff = long_seq_length_cutoff;
    // SOffset heuristics_length_cutoff = std::numeric_limits<SOffset>::max();

    if (max_dist < std::numeric_limits<SOffset>::max()) {
        max_diag_diff = std::numeric_limits<UDiag>::max();
    }

    if (max_diag_diff < std::numeric_limits<UDiag>::max()
        || max_dist < std::numeric_limits<SOffset>::max()) {
        std::cout << "Heuristics"
                  << "\tActivation min. length: " << heuristics_length_cutoff
                  << "\tMax diag diff: " << max_diag_diff << "\tMax dist: " << max_dist
                  << std::endl;
    }

    Offset long_indel_cutoff = 25;

    auto query_rc = reverse_complement(query);

    std::cout << "Aligning: " << theader << " " << target.size() << " -> " << qheader
              << " " << query.size() << "\t" << target_diag << "\n";

    assert(argc > 4);
    int nthreads = std::atoi(argv[4]);
    int nthreads_chaining = nthreads;

    // Feature toggles
    assert(argc > 5);
    bool check_inversions = std::atoi(argv[5]);

    bool check_rc = false;
    if (argc > 6) {
        check_rc = std::atoi(argv[6]);
    }

    bool force_less_heuristics = false;
    if (argc > 7) {
        force_less_heuristics = std::atoi(argv[7]);
    }

    auto check_if_heuristics = [&force_less_heuristics,
                                &heuristics_length_cutoff](SOffset qlen, SOffset tlen) -> bool {
        return !force_less_heuristics && (qlen + tlen > heuristics_length_cutoff);
    };

    assert(argc > 6);
    std::string chain_out(argv[6]);

    std::cout << "Parsing MUMmer\n";
    assert(argc > 3);
    std::ifstream finn(argv[3]);

    auto mummer_ranges = read_mummer(finn, target, theader, query, qheader, query_rc);
    finn.close();

    auto [min_mum_length, total_mum_length] = std::accumulate(
            mummer_ranges.begin(), mummer_ranges.end(),
            std::make_pair<SOffset, SOffset>(std::numeric_limits<SOffset>::max(), 0),
            [](auto&& acc, const auto& range) {
                const auto& [rbegin, rend, rrc, qbegin, qend, qrc, left_trim, right_trim]
                        = range;
                assert(rend >= rbegin);
                assert(rend - rbegin == qend - qbegin);
                SOffset mum_length = rend - rbegin;
                acc.first = std::min(acc.first, mum_length);
                acc.second += mum_length;
                return acc;
            });

    std::cout << "Min. MUM length: " << min_mum_length << "\tAvg. MUM length: "
              << static_cast<double>(total_mum_length) / mummer_ranges.size() << "\n";

    if (min_mum_length < k) {
        throw std::runtime_error("Min. MUM length too short");
    }

    ChainScoreModel chain_score_model(min_mum_length, score_model.gap_switch, exp_mismatch_frac_between_mum_bp, max_gap);
    std::cout << chain_score_model << "\n";

    bool is_rc = false;

    std::cout << "Chaining " << mummer_ranges.size() << " MUMs\n";
    auto [best_chain, chain_score]
            = chain_ranges(target, query, query_rc, mummer_ranges, chain_score_model,
                           check_inversions, nthreads_chaining,
                           true /* show_progress */);

    if (check_rc) {
        auto rc_mummer_ranges = rc_ranges(mummer_ranges, query.size());
        std::cout << "Chaining " << rc_mummer_ranges.size() << " rev-comp MUMs\n";
        auto [best_chain_rc, chain_score_rc]
                = chain_ranges(target, query_rc, query, rc_mummer_ranges, chain_score_model,
                               check_inversions, nthreads_chaining,
                               true /* show_progress */);

        if (best_chain_rc.size() && chain_score_rc > chain_score) {
            std::cout << "Best chain is on rev-comp of query: " << chain_score_rc << " > "
                    << chain_score << "\n";
            is_rc = true;
            std::swap(query_rc, query);
            std::swap(chain_score_rc, chain_score);
            std::swap(best_chain_rc, best_chain);
            std::swap(rc_mummer_ranges, mummer_ranges);
        }
    }

    if (best_chain.empty()) {
        std::cout << "No alignment\n";
        return 0;
    }

    auto [inv_starts, inv_ends] = compute_invs(target, query, query_rc, best_chain);

    reseed_large_gaps(target, query, query_rc,
                      score_model,
                      best_chain,
                      inv_starts, inv_ends,
                      k, check_if_heuristics,
                      max_gap, exp_mismatch_frac_between_mum_bp,
                      check_inversions, nthreads_chaining,
                      true /* show_progress */,
                      false /* show_progress_per_chaining */);

    auto init_storage = [&]() {
        return std::make_tuple(std::vector<std::tuple<int64_t, uint64_t, uint64_t>>(
                                       best_chain.size()),
                               std::vector<Cigar>(best_chain.size()),
                               std::vector<Score>(best_chain.size()),
                               std::vector<int64_t>(best_chain.size()),
                               std::vector<int64_t>(best_chain.size()),
                               std::vector<int64_t>(best_chain.size()),
                               std::vector<int64_t>(best_chain.size()),
                               std::vector<int64_t>(best_chain.size()),
                               std::vector<int64_t>(best_chain.size()),
                               std::vector<std::vector<Ranges>>(best_chain.size()));
    };

    auto [longest_gaps, cigar_parts, scores, identities, matches, inverted, inverted_r, nref, nquery,
          extra_seeds] = init_storage();

    assert(longest_gaps.size() == best_chain.size());
    assert(cigar_parts.size() == best_chain.size());
    assert(scores.size() == best_chain.size());
    assert(identities.size() == best_chain.size());
    assert(matches.size() == best_chain.size());
    assert(inverted.size() == best_chain.size());
    assert(inverted_r.size() == best_chain.size());
    assert(nref.size() == best_chain.size());
    assert(nquery.size() == best_chain.size());
    assert(extra_seeds.size() == best_chain.size());

    assert(best_chain.size() > 2);
    {
        std::ofstream chain_fout(chain_out);
        if (chain_fout) {
            print_chain(chain_fout, theader, target, qheader, query, best_chain, chain_score);
        } else {
            std::cerr << "ERROR: unable to open " << chain_out << " for writing" << std::endl;
            throw std::runtime_error("Failed to open chain output file");
        }
    }

    auto [max_diag, target_abs_diag]
            = compute_diag(target, query, query_rc, best_chain, inv_starts, inv_ends);
    std::cout << "Max diag: " << max_diag << " / " << target_abs_diag << "\t("
              << (target_abs_diag == target_diag ? '+' : '-')
              << ")\n";

    // first pass: handle SNVs and short mismatch regions
    // prefer all mismatches if X * mismatch_p <= 2(gap_open_p + gap_ext_p)
    double larger_gap_cutoff = static_cast<double>(score_model.get_gap_penalty(1) * 2) / score_model.mismatch_p;
    std::cout << "Use Hamming distance for sequences with fewer than " << larger_gap_cutoff << " mismatches\n";
    std::vector<size_t> larger_gaps;
    size_t nrconsumed = 0;
    size_t nqconsumed = 0;
    size_t prealigned = 0;
    for (size_t i = 1; i < best_chain.size(); ++i) {
        Cigar& local_cigar = cigar_parts[i];
        Score& local_score = scores[i];
        bool qorientation_last = best_chain[i - 1].qorientation;
        bool qorientation = best_chain[i].qorientation;

        if (qorientation == qorientation_last) {
            auto [qorientation_last, query_w, target_w, mum_length, qi, ti, exact_match_length]
                    = extract_gap_seqs_continue(target, query, query_rc, best_chain, i);

            if (query_w.size() != target_w.size()) {
                larger_gaps.emplace_back(i);
                continue;
            }

            SOffset qlen = query_w.size();
            SOffset rlen = target_w.size();
            SOffset num = 0;
            SOffset nmismatch = 0;
            char last_op = N_OP;
            for (size_t i = 0; i < query_w.size(); ++i) {
                char cur_op = EQ_OP;
                if (query_w[i] == 'N' || target_w[i] == 'N') {
                    cur_op = MATCH_OP;
                } else if (query_w[i] != target_w[i]) {
                    cur_op = NEQ_OP;
                }
                if (cur_op != last_op) {
                    if (num) {
                        assert(last_op != N_OP);
                        nmismatch += last_op == NEQ_OP ? num : 0;
                        if (nmismatch >= larger_gap_cutoff) {
                            num = 0;
                            break;
                        }
                        local_score += (last_op == NEQ_OP ? score_model.mismatch_s : score_model.match_s) * num;
                        local_cigar.push(last_op, num);
                        num = 0;
                    }
                    last_op = cur_op;
                }
                ++num;
            }
            if (num) {
                assert(last_op != N_OP);
                nmismatch += last_op == NEQ_OP ? num : 0;
                local_score += (last_op == NEQ_OP ? score_model.mismatch_s : score_model.match_s) * num;
                local_cigar.push(last_op, num);
            }
            if (nmismatch >= larger_gap_cutoff) {
                local_score = 0;
                local_cigar.clear();
                larger_gaps.emplace_back(i);
                continue;
            }

            if (exact_match_length) {
                local_score += exact_match_length * score_model.match_s;
                local_cigar.push(EQ_OP, exact_match_length);

                if (qorientation) {
                    local_score += score_model.inv_ext_s * exact_match_length;
                }
            }

            if (qorientation) {
                inverted[i] += exact_match_length + qlen;
                inverted_r[i] += exact_match_length + rlen;
            }

            nref[i] = exact_match_length + rlen;
            nquery[i] = exact_match_length + qlen;
            nrconsumed += nref[i];
            nqconsumed += nquery[i];
            ++prealigned;
        } else if (!qorientation_last) {
            larger_gaps.emplace_back(i);
        }
    }

    std::cout << "Pre-aligned " << prealigned << " / " << best_chain.size() << " regions\n";

    assert(target.size() >= nrconsumed);
    assert(query.size() >= nqconsumed);
    ProgressBar progress_bar(query.size() + target.size() - nrconsumed - nqconsumed, "Global alignment");
    std::mutex mu;

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (size_t j = 0; j < larger_gaps.size(); ++j) {
        size_t i = larger_gaps[j];
        assert(i);

        Cigar& local_cigar = cigar_parts[i];
        Score& local_score = scores[i];
        auto& [gap_length, gap_pos_r, gap_pos_q] = longest_gaps[i];

        bool qorientation = best_chain[i].qorientation;

        auto run_continue = [&](bool qorientation_last, std::string_view query_w,
                                std::string_view target_w, SOffset mum_length, SOffset qi,
                                SOffset ti, SOffset exact_match_length) {
            assert(qorientation == qorientation_last);

            if (query_w.empty() && target_w.empty())
                return;

            size_t skipped = 0;
            while (skipped < query_w.size() && skipped < target_w.size()
                   && (query_w[skipped] == 'N' || target_w[skipped] == 'N')) {
                ++skipped;
            }
            query_w.remove_prefix(skipped);
            target_w.remove_prefix(skipped);
            qi += skipped;
            ti += skipped;
            assert(local_cigar.empty());
            if (skipped) {
                local_cigar.push(MATCH_OP, skipped);
            }

            Score front_eq = std::mismatch(query_w.begin(), query_w.end(), target_w.begin(), target_w.end()).first - query_w.begin();
            if (front_eq) {
                query_w.remove_prefix(front_eq);
                target_w.remove_prefix(front_eq);
                qi += front_eq;
                ti += front_eq;
                local_cigar.push(EQ_OP, front_eq);
            }

            size_t back_skipped = 0;
            while (back_skipped < query_w.size() && back_skipped < target_w.size()
                   && (query_w[query_w.size() - back_skipped - 1] == 'N'
                       || target_w[target_w.size() - back_skipped - 1] == 'N')) {
                ++back_skipped;
            }
            query_w.remove_suffix(back_skipped);
            target_w.remove_suffix(back_skipped);

            local_score += score_model.match_s * static_cast<Score>(skipped + back_skipped + front_eq);

            bool use_heuristics = false;
            bool to_print = false;

            auto print = [&]() {
                SOffset prev_exact_match_length = 0;
                {
                    assert(i);
                    const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation,
                                 left_trim, right_trim] = best_chain[i - 1];
                    assert(rend - rbegin == qend - qbegin);
                    prev_exact_match_length = rend - rbegin;
                }

                std::lock_guard<std::mutex> print_lock(mu);
                std::cout << "Continue run: i: " << i << " / " << best_chain.size() << "\t"
                          << "h: " << use_heuristics << "\t"
                          << "o: " << qorientation_last << "\t"
                          << ti - prev_exact_match_length + 1 << ":"
                          << qi - prev_exact_match_length + 1 << " ("
                          << prev_exact_match_length << ")\t"
                          << front_eq << "\t"
                          << ti + 1 << "-"
                          << ti + target_w.size() << " (" << target_w.size() << ")\t"
                          << qi + 1 << "-" << qi + query_w.size() << " (" << query_w.size()
                          << ")\t" << ti + target_w.size() + exact_match_length << ":"
                          << qi + query_w.size() + exact_match_length << " ("
                          << exact_match_length << ")\t" << mum_length << std::endl;
            };

            if (target_w.empty()) {
                assert(query_w.size());
                local_score += score_model.get_gap_score(query_w.size());
                local_cigar.push(QUERY_CONSUME_OP, query_w.size());
                gap_length = query_w.size();
                if (qorientation) {
                    local_score += score_model.inv_ext_s * static_cast<Score>(query_w.size());
                }
            } else if (query_w.empty()) {
                assert(target_w.size());
                local_score += score_model.get_gap_score(target_w.size());
                local_cigar.push(TARGET_CONSUME_OP, target_w.size());
                gap_length = target_w.size();
            } else {
                use_heuristics = check_if_heuristics(query_w.size(), target_w.size());

                to_print = (use_heuristics
                            || std::max(query_w.size(), target_w.size())
                                    > long_seq_length_cutoff);

                to_print = true;

                if (to_print)
                    print();

                if (has_large_gap(query_w, target_w)) {
                    auto aligner
                            = make_aligner(score_model, use_heuristics ? heuristics_model : main_model);

                    auto [cigar, aln_score] = run_alignment_gaps(*aligner, query_w, target_w, score_model);
                    local_score += aln_score;
                    local_cigar.insert(local_cigar.end(), std::make_move_iterator(cigar.begin()), std::make_move_iterator(cigar.end()));
                } else if (use_heuristics) {
                    Cigar cigar = repeat_aligner(std::string(query_w), std::string(target_w));
                    local_cigar.insert(local_cigar.end(), std::make_move_iterator(cigar.begin()), std::make_move_iterator(cigar.end()));
                    SeqPair view_pair(query_w, target_w);
                    local_score += score_cigar(local_cigar, view_pair, score_model);
                } else {
                    auto aligner = make_aligner(score_model, use_heuristics ? heuristics_model : main_model);
                    if (use_heuristics && max_dist < std::numeric_limits<SOffset>::max()) {
                        aligner->setHeuristicWFadaptive(min_wavefront_length, std::min<SOffset>(max_dist, query_w.size()), 1);
                    }
                    check_lengths(*aligner, score_model.max_pen, query_w.size() + target_w.size());

                    SeqPair view_pair(query_w, target_w);

                    aligner->alignEnd2End(match_char, &view_pair, query_w.size(),
                                         target_w.size());
                    if (aligner->getAlignmentStatus() != wfa::WFAligner::StatusAlgCompleted) {
                        std::cerr << "WARNING: alignment with heuristics failed, rerunning without heuristics\n";
                        aligner->setHeuristicNone();
                        aligner->alignEnd2End(match_char, &view_pair, query_w.size(),
                                         target_w.size());
                    }
                    assert(aligner->getAlignmentStatus() == wfa::WFAligner::StatusAlgCompleted);

                    Score aln_score = 0;
                    Cigar wfa_cigar = aligner->getCIGAR(true);
                    Cigar cigar = cigar_fix_n(wfa_cigar, target_w, query_w);
                    aln_score = score_cigar(cigar, view_pair, score_model);

                    assert(aln_score == score_cigar(cigar, view_pair, score_model));

                    local_score += aln_score;
                    local_cigar.insert(local_cigar.end(), std::make_move_iterator(cigar.begin()), std::make_move_iterator(cigar.end()));
                }
            }

            if (skipped || back_skipped) {
                std::cout << "\n\nShrunk range: " << ti + 1 - skipped << "-"
                          << ti + target_w.size() + back_skipped << " x "
                          << qi + 1 - skipped << "-" << qi + query_w.size() + back_skipped
                          << " to " << ti + 1 << "-" << ti + target_w.size() << " x "
                          << qi + 1 << "-" << qi + query_w.size() << "\n\n"
                          << std::endl;

                if (back_skipped)
                    local_cigar.push(MATCH_OP, back_skipped);
            }

            if (to_print)
                print();

            gap_length = std::max(target_w.size(), query_w.size());
        };

        if (qorientation == best_chain[i - 1].qorientation) {
            auto [qorientation_last, query_w, target_w, mum_length, qi, ti, exact_match_length]
                    = extract_gap_seqs_continue(target, query, query_rc, best_chain, i);

            SOffset rlen = target_w.size();
            SOffset qlen = query_w.size();

            run_continue(qorientation_last, query_w, target_w, mum_length, qi, ti,
                         exact_match_length);

            if (qorientation) {
                inverted[i] = query_w.size();
                inverted_r[i] = target_w.size();
                local_score += score_model.inv_ext_s * inverted[i];
            }

            if (exact_match_length) {
                local_score += exact_match_length * score_model.match_s;
                local_cigar.push(EQ_OP, exact_match_length);

                if (qorientation) {
                    local_score += score_model.inv_ext_s * exact_match_length;
                    inverted[i] += exact_match_length;
                    inverted_r[i] += exact_match_length;
                }
            }


            SOffset r_consumed = rlen + exact_match_length;
            SOffset q_consumed = qlen + exact_match_length;

            nref[i] = r_consumed;
            nquery[i] = q_consumed;

            assert(r_consumed || q_consumed);
            progress_bar += r_consumed + q_consumed;

        } else {
            assert(qorientation);
            assert(i < inv_starts.size());
            assert(i == inv_starts[i]);

            auto [ql_1, query_w_1, query_rc_w_1, target_w_1, mum_length_1, qi_1, ti_1, eml_1]
                = extract_gap_seqs_switch<false>(target, query, query_rc,
                                                 best_chain,
                                                 inv_starts[i] - 1,
                                                 inv_starts[i],
                                                 inv_ends[i],
                                                 inv_ends[i] + 1);
            assert(ql_1 != qorientation);
            assert(target.c_str() + ti_1 == target_w_1.data());
            assert(query.c_str() + qi_1 == query_w_1.data());

            assert(i < inv_ends.size());
            size_t j = inv_ends[i] + 1;
            assert(j < best_chain.size());

            auto [ql_2, query_w_2, query_rc_w_2, target_w_2, mum_length_2, qi_2, ti_2, eml_2]
                = extract_gap_seqs_switch<true>(target, query, query_rc,
                                                best_chain,
                                                inv_starts[i] - 1,
                                                inv_starts[i],
                                                inv_ends[i],
                                                inv_ends[i] + 1);
            assert(ql_2 == qorientation);
            assert(target.c_str() + ti_2 == target_w_2.data());

            assert(query_w_1.size() == query_w_2.size());
            assert(query_rc_w_1.size() == query_rc_w_2.size());

            bool use_heuristics = check_if_heuristics(
                std::max(query_w_1.size(), query_rc_w_1.size()) + std::max(query_w_2.size(), query_rc_w_2.size()),
                target_w_1.size() + target_w_2.size()
            );

            bool to_print = true;

            auto print = [&]() {
                SOffset prev_exact_match_length = 0;
                SOffset ti = 0;
                SOffset qi = 0;
                SOffset next_exact_match_length = 0;
                SOffset rend_next = 0;
                SOffset qend_next = 0;
                {
                    assert(i);
                    const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation,
                                 left_trim, right_trim] = best_chain[i - 1];
                    assert(rend - rbegin == qend - qbegin);
                    prev_exact_match_length = rend - rbegin;
                    ti = rend;
                    qi = qend;
                }
                {
                    assert(j);
                    const auto& [rbegin, rend, rorientation, qbegin, qend, qorientation,
                                 left_trim, right_trim] = best_chain[j];
                    next_exact_match_length = rend - rbegin;
                    rend_next = rend;
                    qend_next = qend;
                }
                std::lock_guard<std::mutex> print_lock(mu);
                std::cout << "Switch run: "
                          << "i: " << i << "\t"
                          << "j: " << j << " / " << best_chain.size() << "\t"
                          << "h: " << use_heuristics << "\t"
                          << ti - prev_exact_match_length + 1 << ":"
                          << qi - prev_exact_match_length + 1 << " ("
                          << prev_exact_match_length << ")\t"
                          << target_w_1.size() << " vs. " << query_w_1.size() << "," << query_rc_w_1.size() << "\t"
                          << eml_1 << "\t"
                          << target_w_2.size() << " vs. " << query_w_2.size() << "," << query_rc_w_2.size() << "\t"
                          << eml_2 << "\t"
                          << rend_next << ":" << qend_next << " ("
                          << next_exact_match_length << ")\t"
                        //   << base_penalty
                          << std::endl;
            };

            if (to_print)
                print();

            auto aligner = make_aligner(score_model, main_model);

            auto &r_consumed_1 = nref[i];
            auto &r_consumed_2 = nref[j];

            Score &score_1 = scores[i];
            Score &score_2 = scores[j];

            Cigar &cigar_1 = cigar_parts[i];
            Cigar &cigar_2 = cigar_parts[j];

            auto &inv_length_1 = inverted[i];
            auto &inv_length_2 = inverted[j];
            auto &inv_length_r_1 = inverted_r[i];
            auto &inv_length_r_2 = inverted_r[j];

            auto &q_consumed_1 = nquery[i];
            auto &q_consumed_2 = nquery[j];

            SOffset inv_len;
            SOffset inv_len_r;

            std::tie(score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
                     score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
                     inv_len, inv_len_r)
                = run_alignment(*aligner, score_model,
                                query, query_rc, target,
                                query_w_1, query_rc_w_1, target_w_1,
                                query_w_2, query_rc_w_2, target_w_2,
                                heuristics_length_cutoff,
                                use_heuristics ? min_wavefront_length : max_diag_width,
                                use_heuristics ? max_dist : max_offset);

            r_consumed_1 += eml_1;
            r_consumed_2 += eml_2;
            q_consumed_1 += eml_1;
            q_consumed_2 += eml_2;

            score_1 += (score_model.match_s + score_model.inv_ext_s) * eml_1;
            score_2 += score_model.match_s * eml_2;

            inv_length_1 += eml_1;
            assert(inv_len >= inv_length_1 + inv_length_2);

            inv_length_r_1 += eml_1;
            assert(inv_len_r >= inv_length_r_1 + inv_length_r_2);

            if (eml_1)
                cigar_1.push(EQ_OP, eml_1);

            if (eml_2)
                cigar_2.push(EQ_OP, eml_2);

            if (to_print)
                print();

            progress_bar += r_consumed_1 + r_consumed_2 + q_consumed_1 + q_consumed_2;
        }
    }

    Score score = std::accumulate(scores.begin(), scores.end(), Score(0ll));

    std::cout << "query is rev-comp: " << is_rc << std::endl;

    std::string query_check;

    Offset r_pos = 0;
    Offset q_pos = 0;
    Offset neq = 0;
    Offset nmatch = 0;
    Offset n_ref = 0;
    Offset n_qry = 0;
    Offset n_longindel_ref = 0;
    Offset n_longindel_qry = 0;

    auto push_op = [&](char last_op, int64_t op_len) {
        assert(q_pos <= query_check.size());
        assert(op_len > 0);
        switch (last_op) {
            case EQ_OP:
            case NEQ_OP:
            case MATCH_OP: {
                assert(r_pos + op_len <= target.size());
                assert(q_pos + op_len <= query.size());

                if (q_pos == query_check.size())
                    query_check += query.substr(q_pos, op_len);

                assert(q_pos + op_len <= query_check.size());

                std::string_view target_w(target.c_str() + r_pos, op_len);
                std::string_view query_w(query_check.c_str() + q_pos, op_len);
                assert(last_op != NEQ_OP || std::equal(target_w.begin(), target_w.end(), query_w.begin(),
                                    [](char a, char b) { return a != b; }));

                if (last_op != MATCH_OP) {
                    assert(target_w.find('N') == std::string_view::npos);
                    assert(query_w.find('N') == std::string_view::npos);
                    nmatch += op_len;
                } else {
                    assert(std::equal(target_w.begin(), target_w.end(), query_w.begin(),
                                      [](char a, char b) { return a == 'N' || b == 'N'; }));
                    n_ref += std::count(target_w.begin(), target_w.end(), 'N');
                    n_qry += std::count(query_w.begin(), query_w.end(), 'N');
                }

                if (last_op == EQ_OP) {
                    assert(target_w == query_w);
                    neq += op_len;
                }

                r_pos += op_len;
                q_pos += op_len;
            } break;
            case TARGET_CONSUME_OP: {
                assert(r_pos + op_len <= target.size());

                std::string_view target_w(target.c_str() + r_pos, op_len);
                Offset n = std::count(target_w.begin(), target_w.end(), 'N');
                n_ref += n;

                if (op_len - n > long_indel_cutoff) {
                    n_longindel_ref += op_len - n;
                }
                r_pos += op_len;
            } break;
            case QUERY_CONSUME_OP: {
                assert(q_pos + op_len <= query.size());

                if (q_pos == query_check.size())
                    query_check += query.substr(q_pos, op_len);

                assert(q_pos + op_len <= query_check.size());

                std::string_view query_w(query_check.c_str() + q_pos, op_len);
                Offset n = std::count(query_w.begin(), query_w.end(), 'N');
                n_qry += n;
                if (op_len - n > long_indel_cutoff) {
                    n_longindel_qry += op_len - n;
                }
                q_pos += op_len;
            } break;
            case INV_OP: {
                assert(check_inversions);
                assert(q_pos + op_len <= query.size());
                query_check += query_rc.substr(query_rc.size() - q_pos - op_len, op_len);
                assert(is_reverse_complement(std::string_view(query_check.c_str() + q_pos, op_len),
                                             std::string_view(query.c_str() + q_pos, op_len)));
            }
        }

        std::cout << op_len << last_op;
    };

    for (size_t i = 0; i < cigar_parts.size(); ++i) {
        cigar_caller(cigar_parts[i], [&](char c, Offset num, Offset /* rel_r_pos */, Offset /* rel_q_pos */) {
            push_op(c, num);
        }, nref[i]);
    }

    std::cout << std::endl;

    assert(query_check.size() == query.size());

    std::cout << neq << " " << nmatch << " " << n_ref << " " << n_qry << " " << n_longindel_ref << " " << n_longindel_qry << "\n";

    Offset rnon = target.size() - n_ref;
    Offset qnon = query.size() - n_qry;
    Offset min_non = std::min(rnon, qnon);
    Offset max_non = std::max(rnon, qnon);

    Offset rnon_nolong = rnon - n_longindel_ref;
    Offset qnon_nolong = qnon - n_longindel_qry;
    Offset min_nolong = std::min(rnon_nolong, qnon_nolong);
    Offset max_nolong = std::max(rnon_nolong, qnon_nolong);
    std::cout << static_cast<double>(neq*100) / max_non << "-"
              << static_cast<double>(neq*100) / min_non << "\t"
              << static_cast<double>(nmatch*100) / max_non << "-"
              << static_cast<double>(nmatch*100) / min_non << "\t"
              << static_cast<double>(neq*100) / nmatch << "\t"
              << static_cast<double>(neq*100) / max_nolong << "-"
              << static_cast<double>(neq*100) / min_nolong
              << "\n";

    assert(r_pos == target.size());
    assert(q_pos == query.size());

    return 0;
}
