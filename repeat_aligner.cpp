#include "repeat_aligner.hpp"

#include "cigar.hpp"

#include <mummer/sparseSA.hpp>
#include <tandem_aligner.hpp>

std::pair<Score, std::string>
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
            cigar = cigar_fix_n(base_cigar, target, query);
        }
    } else if (query.size()) {
        cigar = std::to_string(query.size()) + "I";
    } else if (target.size()) {
        cigar = std::to_string(target.size()) + "D";
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

    return std::make_pair(score, std::move(cigar));
}

std::string repeat_aligner(const std::string &query,
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

    ta.AssignMismatches(
        ta_cigar,
        target,
        query
    );

    std::ostringstream sout;
    sout << ta_cigar;
    return cigar_fix_n(sout.str(), target, query);
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
        std::transform(str.begin(), str.end(), str.begin(), [nucleotides_only](char c) {
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

    if (query.size() >= min_length) {
        sa.findMAM_each(query_str.data(), query_str.size(), min_length, false, [&](const auto& m) {
            assert(m.len);
            Ranges mum(m.ref + tbegin, m.ref + m.len + tbegin, false,
                       m.query + qbegin, m.query + m.len + qbegin, false);
            assert(mum.check_equal(target, query, query_rc));
            callback(std::move(mum));
        });
    }

    if (query_rc.size() >= min_length) {
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
