#include "repeat_aligner.hpp"

#include "cigar.hpp"

#include <mummer/sparseSA.hpp>
#include <tandem_aligner.hpp>

std::pair<Score, Cigar>
get_alignment(wfa::WFAligner& aligner,
                const ScoreModel &score_model,
                std::string_view query,
                std::string_view target,
                bool penalty_to_score,
                SOffset heuristics_length_cutoff,
                Diag min_k,
                Diag max_k) {
    SeqPair view_pair(query, target);
    Score n_penalty = INT32_MIN;
    Cigar cigar;
    if (query.size() && target.size()) {
        if (static_cast<SOffset>(query.size() + target.size()) >= heuristics_length_cutoff) {
            cigar = repeat_aligner(aligner, std::string(query), std::string(target));
        } else {
            if (min_k == min_diag && max_k == max_diag) {
                aligner.setHeuristicNone();
            } else {
                aligner.setHeuristicBandedStatic(-max_k - 1, -min_k + 1);
            }
            aligner.alignEnd2End(match_char, &view_pair, query.size(), target.size());
            if (aligner.getAlignmentStatus() != wfa::WFAligner::StatusAlgCompleted) {
                assert(min_k != min_diag || max_k != max_diag);
                std::cerr << "WARNING: rerunning alignment without heuristics\n";
                aligner.setHeuristicNone();
                aligner.alignEnd2End(match_char, &view_pair, query.size(),
                                    target.size());
            }
            assert(aligner.getAlignmentStatus() == wfa::WFAligner::StatusAlgCompleted);
            n_penalty = aligner.getAlignmentScore();
            std::string base_cigar = aligner.getCIGAR(true);
            cigar = cigar_fix_n(Cigar(base_cigar), target, query);
        }
    } else if (query.size()) {
        cigar = Cigar(QUERY_CONSUME_OP, query.size());
    } else if (target.size()) {
        cigar = Cigar(TARGET_CONSUME_OP, target.size());
    } else {
        n_penalty = 0;
    }

    assert((query.empty() && target.empty()) || !cigar.empty());
    assert(cigar == cigar_fix_n(cigar, target, query));
    Score score;
    if (n_penalty != INT32_MIN) {
        score = penalty_to_score
            ? score_model.penalty_to_score(-n_penalty, query.size(), target.size())
            : -n_penalty;
    } else {
        score = score_cigar(cigar, view_pair, score_model, !penalty_to_score);
    }

    return std::make_pair(score, std::move(cigar));
}

Cigar repeat_aligner(wfa::WFAligner& aligner,
                     const std::string &query,
                     const std::string &target) {
    tandem_aligner::Cigar ta_cigar;
    std::queue<tandem_aligner::MinSeqTask> queue;
    queue.push({
        ta_cigar.begin(),
        0, (int64_t) target.size(),
        0, (int64_t) query.size()
    });
    int max_freq = 50;
    logging::Logger logger;
    tandem_aligner::TandemAligner ta(logger,
        "./tmp",
        max_freq,
        false,
        false,
        false);
    while (queue.size()) {
        ta.RunTask(queue, ta_cigar,
                   target,
                   query,
                   false, false);
        queue.pop();
    }

    ta_cigar.AssertValidity(target, query);
    // ta.AssignMismatches(
    //     ta_cigar,
    //     target,
    //     query
    // );
    // std::ostringstream sout;
    // sout << ta_cigar;
    // return cigar_fix_n(Cigar(sout.str()), target, query);

    Cigar cigar;

    auto push_meqx = [&](std::string_view target_w, std::string_view query_w) {
        assert(target_w.size() == query_w.size());
        for (size_t i = 0; i < target_w.size(); ++i) {
            if (target_w[i] == 'N' || query_w[i] == 'N') {
                cigar.push(MATCH_OP, 1, true);
            } else {
                cigar.push(target_w[i] == query_w[i] ? EQ_OP : NEQ_OP, 1, true);
            }
        }
    };

    if (ta_cigar.Size() >= 2) {
        size_t target_i = 0;
        size_t query_i = 0;
        auto last = ta_cigar.begin();
        auto cur = ta_cigar.begin();
        for (++cur; cur != ta_cigar.end(); ++cur, ++last) {
            size_t length = last->length;
            switch (last->mode) {
                case tandem_aligner::CigarMode::M:
                case tandem_aligner::CigarMode::X: {
                    assert(target_i + length <= target.size());
                    assert(query_i + length <= query.size());
                    std::string_view target_w(target.c_str() + target_i, length);
                    std::string_view query_w(query.c_str() + query_i, length);
                    push_meqx(target_w, query_w);
                    target_i += length;
                    query_i += length;
                } break;
                case tandem_aligner::CigarMode::I:
                case tandem_aligner::CigarMode::D: {
                    switch (cur->mode) {
                        case tandem_aligner::CigarMode::I:
                        case tandem_aligner::CigarMode::D: {
                            // replace both with an alignment
                            size_t query_delta = (last->mode == tandem_aligner::CigarMode::I) * last->length
                                                    + (cur->mode == tandem_aligner::CigarMode::I) * cur->length;
                            size_t target_delta = (last->mode == tandem_aligner::CigarMode::D) * last->length
                                                    + (cur->mode == tandem_aligner::CigarMode::D) * cur->length;
                            assert(query_i + query_delta <= query.size());
                            assert(target_i + target_delta <= target.size());
                            std::string_view target_w(target.c_str() + target_i, target_delta);
                            std::string_view query_w(query.c_str() + query_i, query_delta);
                            SeqPair view_pair(query_w, target_w);
                            aligner.setHeuristicNone();
                            aligner.alignEnd2End(match_char, &view_pair, query_w.size(), target_w.size());
                            assert(aligner.getAlignmentStatus() == wfa::WFAligner::StatusAlgCompleted);
                            std::string base_cigar = aligner.getCIGAR(true);
                            auto local_cigar = cigar_fix_n(Cigar(base_cigar), target_w, query_w);
                            cigar.insert(
                                cigar.end(),
                                std::make_move_iterator(local_cigar.begin()),
                                std::make_move_iterator(local_cigar.end())
                            );
                            query_i += query_delta;
                            target_i += target_delta;
                            ++last;
                            ++cur;
                        } break;
                        default: {
                            if (last->mode == tandem_aligner::CigarMode::I) {
                                assert(query_i + length <= query.size());
                                cigar.push(QUERY_CONSUME_OP, length, true);
                                query_i += length;
                            } else if (last->mode == tandem_aligner::CigarMode::D) {
                                assert(target_i + length <= target.size());
                                cigar.push(TARGET_CONSUME_OP, length, true);
                                target_i += length;
                            } else {
                                assert(false && "This should not happen");
                            }
                        } break;
                    }
                } break;
            }

            if (cur == ta_cigar.end()) {
                ++last;
                assert(last == ta_cigar.end());
                break;
            }
        }

        if (last != ta_cigar.end()) {
            size_t length = last->length;
            switch (last->mode) {
                case tandem_aligner::CigarMode::M:
                case tandem_aligner::CigarMode::X: {
                    assert(target_i + length <= target.size());
                    assert(query_i + length <= query.size());
                    std::string_view target_w(target.c_str() + target_i, length);
                    std::string_view query_w(query.c_str() + query_i, length);
                    push_meqx(target_w, query_w);
                    target_i += length;
                    query_i += length;
                } break;
                case tandem_aligner::CigarMode::I: {
                    assert(query_i + length <= query.size());
                    cigar.push(QUERY_CONSUME_OP, length, true);
                    query_i += length;
                } break;
                case tandem_aligner::CigarMode::D: {
                    assert(target_i + length <= target.size());
                    cigar.push(TARGET_CONSUME_OP, length, true);
                    target_i += length;
                } break;
            }
        }

    } else if (ta_cigar.Size()) {
        size_t length = ta_cigar.begin()->length;
        switch (ta_cigar.begin()->mode) {
            case tandem_aligner::CigarMode::M:
            case tandem_aligner::CigarMode::X: {
                assert(target.size() == length);
                assert(query.size() == length);
                push_meqx(target, query);
            } break;
            case tandem_aligner::CigarMode::I: {
                assert(query.size() == length);
                cigar.push(QUERY_CONSUME_OP, length, true);
            } break;
            case tandem_aligner::CigarMode::D: {
                assert(target.size() == length);
                cigar.push(TARGET_CONSUME_OP, length, true);
            } break;
        }
    }

    return cigar;
}

void call_mums(std::string_view query,
               std::string_view query_rc,
               std::string_view target,
               size_t min_length,
               const std::function<void(Ranges&&)> &callback,
               SOffset qbegin,
               SOffset qend,
               SOffset qrcbegin,
               SOffset qrcend,
               SOffset tbegin,
               SOffset tend) {
    if (target.size() < min_length)
        return;

    if (qend == -1)
        qend = query.size();

    if (qrcend == -1)
        qrcend = query_rc.size();

    if (tend == -1)
        tend = target.size();

    std::vector<std::string> descr;
    std::vector<long> startpos;

    descr.emplace_back();
    startpos.emplace_back(0);

    auto prepare_str = [](std::string_view strview, bool nucleotides_only = true) {
        std::string str(strview);
        std::transform(str.begin(), str.end(), str.begin(), [nucleotides_only](unsigned char c) {
            c = std::tolower(c);
            if (nucleotides_only) {
                switch (c) {
                    case 'a':
                    case 't':
                    case 'g':
                    case 'c': break;
                    default: c = '~';
                }
            }
            return c;
        });
        return str;
    };

    std::string target_str = prepare_str(target.substr(tbegin, tend - tbegin), false);
    std::string query_str = prepare_str(query.substr(qbegin, qend - qbegin));
    std::string query_rc_str = prepare_str(query_rc.substr(qrcbegin, qrcend - qrcbegin));

    // Build the enhanced sparse suffix array over the reference. create_auto
    // derives the sparseness (K, sparseMult, child/suffix-link tables) from
    // min_len. NOTE: the object stores a raw pointer into ref's buffer, so ref
    // must outlive `sa`.
    long K = 1;
    bool suflink = (K < 4);
    bool child = (K >= 4);
    int sparseMult = 1;

    if (suflink && !child) {
        sparseMult = 1;
    } else if (K >= 4) {
        sparseMult = std::max<int>((static_cast<int>(min_length) - 10) / K, 1);
    } else {
        sparseMult = std::max<int>((static_cast<int>(min_length) - 12) / K, 1);
    }

    int kmer = std::max<int>(0, std::min<int>(10, static_cast<int>(min_length) - sparseMult * K + 1));

    bool _4column = false;
    bool printSubstring = false;
    bool nucleotidesOnly = true;

    mummer::mummer::sparseSAMatch sa(target_str, descr, startpos,
                                     _4column, K,
                                     suflink, child, kmer > 0,
                                     sparseMult, kmer, printSubstring, nucleotidesOnly);
    sa.construct();

    if (query_str.size() >= min_length) {
        sa.findMAM_each(query_str.data(), query_str.size(), min_length, false, [&](const auto& m) {
            assert(m.len);
            Ranges mum(m.ref + tbegin, m.ref + m.len + tbegin, false,
                       m.query + qbegin, m.query + m.len + qbegin, false);
            assert(mum.check_equal(target, query, query_rc));
            callback(std::move(mum));
        });
    }

    if (query_rc_str.size() >= min_length) {
        sa.findMAM_each(query_rc_str.data(), query_rc_str.size(), min_length, false, [&](const auto& m) {
            assert(m.len);
            Ranges rc_mum(m.ref + tbegin, m.ref + m.len + tbegin, false,
                          query_rc.size() - (m.query + m.len + qrcbegin),
                          query_rc.size() - (m.query + qrcbegin),
                          true);
            assert(rc_mum.check_equal(target, query, query_rc));
            callback(std::move(rc_mum));
        });
    }
}
