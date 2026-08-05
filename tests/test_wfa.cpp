#include <gtest/gtest.h>

#include "helpers.hpp"
#include "chaining.hpp"
#include "repeat_aligner.hpp"
#include "wfa_switch.hpp"

#include <string>
#include <string_view>

TEST(WFATest, CheckAlignment) {
    std::string target = "ACGTAC";
    std::string query  = "ACGTACGTAC";
    SeqPair view_pair(query, target);

    ASSERT_GT(query.size(), target.size());
    Cigar cigar_truth(EQ_OP, target.size());
    cigar_truth.push(QUERY_CONSUME_OP, query.size() - target.size());

    ScoreModel score_model(1, -9, -16, -2, -41, -1, -41, 0);
    auto aligner = make_aligner(score_model);

    aligner->alignEnd2End(query, target);
    ASSERT_EQ(wfa::WFAligner::StatusAlgCompleted, aligner->getAlignmentStatus());
    EXPECT_EQ(cigar_truth, aligner->getCIGAR(true));

    aligner->alignEnd2End(match_char, &view_pair, query.size(), target.size());
    ASSERT_EQ(wfa::WFAligner::StatusAlgCompleted, aligner->getAlignmentStatus());
    EXPECT_EQ(cigar_truth, aligner->getCIGAR(true));

    auto [score, cigar] = get_alignment(*aligner, score_model, query, target);
    EXPECT_EQ(cigar_truth, cigar);
    EXPECT_EQ(score, score_cigar(cigar, view_pair, score_model));
}

TEST(WFATest, ZeroLength) {
    std::string target;
    std::string query;
    std::string_view target_w(target);
    std::string_view query_w(query);

    ScoreModel score_model;

    WFAIterator<> wfa_it(score_model,
                         query_w.begin(), query_w.end(),
                         target_w.begin(), target_w.end());
    wfa_it.extend();

    for (Penalty p = 1; p <= score_model.max_pen + 1; ++p) {
        EXPECT_EQ(0u, wfa_it.next());
        EXPECT_EQ(p, wfa_it.get_p());
    }

    EXPECT_TRUE(wfa_it.empty());
}

TEST(WFATest, DualZeroLength) {
    std::string target = "TTAACC";
    std::string query = "TTTTCC";
    std::string query_rc = reverse_complement(query);

    std::string_view target_1(target.c_str() + 2, 0);
    std::string_view query_1(query.c_str() + 2, 0);
    std::string_view query_rc_1(query_rc.c_str() + 2, 0);
    std::string_view target_2(target.c_str() + 4, 0);
    std::string_view query_2(query_rc.c_str() + 4, 0);
    std::string_view query_rc_2(query.c_str() + 4, 0);

    ScoreModel score_model;

    auto aligner = make_aligner(score_model);

    auto [score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
          score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
          inv_len, inv_len_r]
        = run_alignment(*aligner, score_model,
                        query, query_rc, target,
                        query_1, query_rc_1, target_1,
                        query_2, query_rc_2, target_2);

    EXPECT_EQ(0, score_1);
    EXPECT_EQ(Cigar(INV_OP, 2), cigar_1);
    EXPECT_EQ(0, r_consumed_1);
    EXPECT_EQ(0, q_consumed_1);
    EXPECT_EQ(0, inv_length_1);
    EXPECT_EQ(0, inv_length_r_1);

    EXPECT_EQ(0, score_2);
    EXPECT_EQ(Cigar(), cigar_2);
    EXPECT_EQ(0, r_consumed_2);
    EXPECT_EQ(0, q_consumed_2);
    EXPECT_EQ(0, inv_length_2);
    EXPECT_EQ(0, inv_length_r_2);

    EXPECT_EQ(2, inv_len);
    EXPECT_EQ(2, inv_len_r);
}

TEST(WFATest, SharedMiddle) {
    std::string shared_prefix = "TTTAGGATTGCAGATGATTTTGCTTTTAGTTATATTTTCCACACACTGTCCAGAATTACC";
    std::string target_middle = "AGATGTTATACACAAAATGTCGAGCAGACTTCACTTGAGAAAACATTCTTTCTCAGGATT";
    std::string shared_suffix = "CCAGTCACAACCTGCAATTGTGCGGCAAAGTGCAGCAAAGTATTTCCTAAATCTGCAAAC";
    std::string query_middle = reverse_complement(target_middle);

    std::string target = shared_prefix + target_middle + shared_suffix;
    std::string query = shared_prefix + query_middle + shared_suffix;
    std::string query_rc = reverse_complement(query);

    ASSERT_EQ(target.size(), query.size());

    std::vector<Ranges> ranges;
    ranges.emplace_back(0, 0, false, 0, 0, false, 0, 0);
    ranges.emplace_back(0, shared_prefix.size() + 1, false,
                        0, shared_prefix.size() + 1, false,
                        0, 0);
    ranges.emplace_back(shared_prefix.size(), shared_prefix.size() + target_middle.size() + 1, false,
                        shared_prefix.size(), shared_prefix.size() + target_middle.size() + 1, true,
                        0, 0);
    ranges.emplace_back(shared_prefix.size() + target_middle.size() - 1, target.size(), false,
                        shared_prefix.size() + target_middle.size() - 1, query.size(), false,
                        0, 0);
    ranges.emplace_back(target.size(), target.size(), false, query.size(), query.size(), false, 0, 0);

    auto [ql_1, query_w_1, query_rc_w_1, target_w_1, mum_length_1, qi_1, ti_1, eml_1]
        = extract_gap_seqs_switch<false>(target, query, query_rc, ranges, 1, 2, 2, 3);
    ASSERT_FALSE(ql_1);

    auto [ql_2, query_w_2, query_rc_w_2, target_w_2, mum_length_2, qi_2, ti_2, eml_2]
        = extract_gap_seqs_switch<true>(target, query, query_rc, ranges, 1, 2, 2, 3);
    ASSERT_TRUE(ql_2);

    ASSERT_EQ(query_w_1.size(), query_w_2.size());
    ASSERT_EQ(query_rc_w_1.size(), query_rc_w_2.size());

    ScoreModel score_model;

    auto aligner = make_aligner(score_model);

    auto [score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
          score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
          inv_len, inv_len_r]
        = run_alignment(*aligner, score_model,
                        query, query_rc, target,
                        query_w_1, query_rc_w_1, target_w_1,
                        query_w_2, query_rc_w_2, target_w_2);

    EXPECT_EQ(0, score_1);
    EXPECT_EQ(Cigar(INV_OP, target_middle.size() - 2), cigar_1);
    EXPECT_EQ(0, r_consumed_1);
    EXPECT_EQ(0, q_consumed_1);
    EXPECT_EQ(0, inv_length_1);
    EXPECT_EQ(0, inv_length_r_1);

    EXPECT_EQ(0, score_2);
    EXPECT_EQ(Cigar(), cigar_2);
    EXPECT_EQ(0, r_consumed_2);
    EXPECT_EQ(0, q_consumed_2);
    EXPECT_EQ(0, inv_length_2);
    EXPECT_EQ(0, inv_length_r_2);

    EXPECT_EQ(target_middle.size() - 2, inv_len);
    EXPECT_EQ(target_middle.size() - 2, inv_len_r);
}

TEST(WFATest, EmptyAltQueriesEq) {
    std::string shared_prefix = "TT";
    std::string shared_suffix = "AA";
    std::string shared_middle = "CC";

    std::string t1 = "ATGC";
    std::string t2 = reverse_complement(t1);
    std::string q1 = "ATGC";
    std::string qrc1;
    std::string q2 = reverse_complement(q1);
    std::string qrc2 = reverse_complement(qrc1);

    ASSERT_TRUE(is_reverse_complement(q1, q2));
    ASSERT_TRUE(is_reverse_complement(qrc1, qrc2));

    ASSERT_EQ(t1, q1);

    std::string target = shared_prefix + t1 + shared_middle + t2 + shared_suffix;
    std::string query = shared_prefix + q1 + reverse_complement(shared_middle) + qrc2 + shared_suffix;
    std::string query_rc = reverse_complement(query);

    std::string_view target_1(target.c_str() + shared_prefix.size(), t1.size());
    ASSERT_EQ(t1, target_1);

    std::string_view query_1(query.c_str() + shared_prefix.size(), q1.size());
    ASSERT_EQ(q1, query_1);

    std::string_view query_rc_1(query_rc.c_str() + shared_suffix.size(), qrc1.size());
    std::string_view target_2(target_1.data() + t1.size() + shared_middle.size(), t2.size());
    std::string_view query_2(query_rc_1.data() + qrc1.size() + shared_middle.size(), q2.size());
    std::string_view query_rc_2(query_1.data() + q1.size() + shared_middle.size(), qrc2.size());

    ASSERT_EQ(qrc1, query_rc_1);
    ASSERT_EQ(t2, target_2);
    ASSERT_EQ(q2, query_2);
    ASSERT_EQ(qrc2, query_rc_2);

    Cigar cigar_1_truth(EQ_OP, target_1.size());
    cigar_1_truth.push(INV_OP, shared_middle.size());

    Cigar cigar_2_truth(TARGET_CONSUME_OP, target_2.size());

    ScoreModel score_model;

    auto aligner = make_aligner(score_model);

    auto [score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
          score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
          inv_len, inv_len_r]
        = run_alignment(*aligner, score_model,
                        query, query_rc, target,
                        query_1, query_rc_1, target_1,
                        query_2, query_rc_2, target_2);

    EXPECT_EQ(score_model.match_s * target_1.size(), score_1);
    EXPECT_EQ(cigar_1_truth, cigar_1);
    EXPECT_EQ(target_1.size(), r_consumed_1);
    EXPECT_EQ(query_1.size(), q_consumed_1);
    EXPECT_EQ(0, inv_length_1);
    EXPECT_EQ(0, inv_length_r_1);

    EXPECT_EQ(score_model.get_gap_score(target_2.size()), score_2);
    EXPECT_EQ(cigar_2_truth, cigar_2);
    EXPECT_EQ(target_2.size(), r_consumed_2);
    EXPECT_EQ(0, q_consumed_2);
    EXPECT_EQ(0, inv_length_2);
    EXPECT_EQ(0, inv_length_r_2);

    EXPECT_EQ(shared_middle.size(), inv_len);
    EXPECT_EQ(shared_middle.size(), inv_len_r);
}

TEST(WFATest, EmptyAltQueriesNeq) {
    std::string shared_prefix = "TT";
    std::string shared_suffix = "AA";
    std::string shared_middle = "CC";

    std::string t1 = "ATGC";
    std::string q1 = "ATGC";
    std::string qrc1;
    std::string t2 = "ATGC";
    std::string q2 = reverse_complement(q1);
    std::string qrc2 = reverse_complement(qrc1);

    ASSERT_TRUE(is_reverse_complement(q1, q2));
    ASSERT_TRUE(is_reverse_complement(qrc1, qrc2));

    std::string target = shared_prefix + t1 + shared_middle + t2 + shared_suffix;
    std::string query = shared_prefix + q1 + reverse_complement(shared_middle) + qrc2 + shared_suffix;
    std::string query_rc = reverse_complement(query);

    std::string_view target_1(target.c_str() + shared_prefix.size(), t1.size());
    ASSERT_EQ(t1, target_1);

    std::string_view query_1(query.c_str() + shared_prefix.size(), q1.size());
    ASSERT_EQ(q1, query_1);

    std::string_view query_rc_1(query_rc.c_str() + shared_suffix.size(), qrc1.size());
    std::string_view target_2(target_1.data() + t1.size() + shared_middle.size(), t2.size());
    std::string_view query_2(query_rc_1.data() + qrc1.size() + shared_middle.size(), q2.size());
    std::string_view query_rc_2(query_1.data() + q1.size() + shared_middle.size(), qrc2.size());

    ASSERT_EQ(qrc1, query_rc_1);
    ASSERT_EQ(t2, target_2);
    ASSERT_EQ(q2, query_2);
    ASSERT_EQ(qrc2, query_rc_2);

    Cigar cigar_1_truth(EQ_OP, target_1.size());
    cigar_1_truth.push(INV_OP, shared_middle.size());

    Cigar cigar_2_truth(TARGET_CONSUME_OP, target_2.size());

    ScoreModel score_model;

    auto aligner = make_aligner(score_model);

    auto [score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
          score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
          inv_len, inv_len_r]
        = run_alignment(*aligner, score_model,
                        query, query_rc, target,
                        query_1, query_rc_1, target_1,
                        query_2, query_rc_2, target_2);

    EXPECT_EQ(score_model.match_s * target_1.size(), score_1);
    EXPECT_EQ(cigar_1_truth, cigar_1);
    EXPECT_EQ(target_1.size(), r_consumed_1);
    EXPECT_EQ(query_1.size(), q_consumed_1);
    EXPECT_EQ(0, inv_length_1);
    EXPECT_EQ(0, inv_length_r_1);

    EXPECT_EQ(score_model.get_gap_score(target_2.size()), score_2);
    EXPECT_EQ(cigar_2_truth, cigar_2);
    EXPECT_EQ(target_2.size(), r_consumed_2);
    EXPECT_EQ(0, q_consumed_2);
    EXPECT_EQ(0, inv_length_2);
    EXPECT_EQ(0, inv_length_r_2);

    EXPECT_EQ(shared_middle.size(), inv_len);
    EXPECT_EQ(shared_middle.size(), inv_len_r);
}

TEST(WFATest, EmptyQueries) {
    std::string shared_prefix = "TT";
    std::string shared_suffix = "AA";
    std::string shared_middle = "CC";

    std::string t1 = "ATGC";
    std::string q1;
    std::string qrc1;
    std::string t2 = "ATGC";
    std::string q2;
    std::string qrc2;

    ASSERT_TRUE(is_reverse_complement(q1, q2));
    ASSERT_TRUE(is_reverse_complement(qrc1, qrc2));

    std::string target = shared_prefix + t1 + shared_middle + t2 + shared_suffix;
    std::string query = shared_prefix + q1 + reverse_complement(shared_middle) + qrc2 + shared_suffix;
    std::string query_rc = reverse_complement(query);

    std::string_view target_1(target.c_str() + shared_prefix.size(), t1.size());
    ASSERT_EQ(t1, target_1);

    std::string_view query_1(query.c_str() + shared_prefix.size(), q1.size());
    ASSERT_EQ(q1, query_1);

    std::string_view query_rc_1(query_rc.c_str() + shared_suffix.size(), qrc1.size());
    std::string_view target_2(target_1.data() + t1.size() + shared_middle.size(), t2.size());
    std::string_view query_2(query_rc_1.data() + qrc1.size() + shared_middle.size(), q2.size());
    std::string_view query_rc_2(query_1.data() + q1.size() + shared_middle.size(), qrc2.size());

    ASSERT_EQ(qrc1, query_rc_1);
    ASSERT_EQ(t2, target_2);
    ASSERT_EQ(q2, query_2);
    ASSERT_EQ(qrc2, query_rc_2);

    Cigar cigar_1_truth(TARGET_CONSUME_OP, target_1.size());
    cigar_1_truth.push(INV_OP, shared_middle.size());

    Cigar cigar_2_truth(TARGET_CONSUME_OP, target_2.size());

    ScoreModel score_model;

    auto aligner = make_aligner(score_model);

    auto [score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
          score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
          inv_len, inv_len_r]
        = run_alignment(*aligner, score_model,
                        query, query_rc, target,
                        query_1, query_rc_1, target_1,
                        query_2, query_rc_2, target_2);

    EXPECT_EQ(score_model.get_gap_score(target_1.size()) + score_model.inv_open_s, score_1);
    EXPECT_EQ(cigar_1_truth, cigar_1);
    EXPECT_EQ(target_1.size(), r_consumed_1);
    EXPECT_EQ(0, q_consumed_1);
    EXPECT_EQ(0, inv_length_1);
    EXPECT_EQ(0, inv_length_r_1);

    EXPECT_EQ(score_model.get_gap_score(target_2.size()), score_2);
    EXPECT_EQ(cigar_2_truth, cigar_2);
    EXPECT_EQ(target_2.size(), r_consumed_2);
    EXPECT_EQ(0, q_consumed_2);
    EXPECT_EQ(0, inv_length_2);
    EXPECT_EQ(0, inv_length_r_2);

    EXPECT_EQ(shared_middle.size(), inv_len);
    EXPECT_EQ(shared_middle.size(), inv_len_r);
}

TEST(WFATest, BothRight) {
    std::string shared_prefix = "TT";
    std::string shared_suffix = "AA";
    std::string shared_middle = "CC";

    // target_1 aligns to nothing, target_2 aligns to query_2 + query_rc_2
    std::string t1 = "CCC";
    std::string t2 = "ACCTTTCT";
    std::string q2 = "ACCT";
    std::string qrc2 = "TTCT";
    std::string q1 = "AGGT";
    std::string qrc1 = "AGAA";

    ASSERT_TRUE(is_reverse_complement(q1, q2));
    ASSERT_TRUE(is_reverse_complement(qrc1, qrc2));

    std::string target = shared_prefix + t1 + shared_middle + t2 + shared_suffix;
    std::string query = shared_prefix + q1 + reverse_complement(shared_middle) + qrc2 + shared_suffix;
    std::string query_rc = reverse_complement(query);

    std::string_view target_1(target.c_str() + shared_prefix.size(), t1.size());
    ASSERT_EQ(t1, target_1);

    std::string_view query_1(query.c_str() + shared_prefix.size(), q1.size());
    ASSERT_EQ(q1, query_1);

    std::string_view query_rc_1(query_rc.c_str() + shared_suffix.size(), qrc1.size());
    std::string_view target_2(target_1.data() + t1.size() + shared_middle.size(), t2.size());
    std::string_view query_2(query_rc_1.data() + qrc1.size() + shared_middle.size(), q2.size());
    std::string_view query_rc_2(query_1.data() + q1.size() + shared_middle.size(), qrc2.size());

    ASSERT_EQ(qrc1, query_rc_1);
    ASSERT_EQ(t2, target_2);
    ASSERT_EQ(q2, query_2);
    ASSERT_EQ(qrc2, query_rc_2);

    ScoreModel score_model;

    auto aligner = make_aligner(score_model);

    auto [score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
          score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
          inv_len, inv_len_r]
        = run_alignment(*aligner, score_model,
                        query, query_rc, target,
                        query_1, query_rc_1, target_1,
                        query_2, query_rc_2, target_2);

    EXPECT_EQ(Cigar("3D6i"), cigar_1);
    EXPECT_EQ(Cigar("4=4="), cigar_2);

    EXPECT_EQ(6, inv_len);
    EXPECT_EQ(6, inv_len_r);
}

TEST(WFATest, BothRightFull) {
    std::string target_1_truth = "CCC";
    std::string target_2_truth = "ACCTTTCT";
    std::string query_2_truth = "";
    std::string query_rc_2_truth = "ACCTTTCT";
    std::string query_1_truth = reverse_complement(query_2_truth);
    std::string query_rc_1_truth = reverse_complement(query_rc_2_truth);

    std::string fw_seed_1 = "ATTGGTT";
    std::string fw_seed_2 = "GTGAGCGATGAC";
    std::string rc_seed_target = "GCATGACGAGCGATC";
    std::string rc_seed_query = reverse_complement(rc_seed_target);

    // target_1 aligns to nothing, target_2 aligns to query_2 + query_rc_2
    std::string target = fw_seed_1 + target_1_truth + rc_seed_target + target_2_truth + fw_seed_2;
    std::string query = fw_seed_1 + rc_seed_query + query_1_truth + query_rc_2_truth + fw_seed_2;
    std::string query_rc = reverse_complement(query);

    std::vector<Ranges> ranges;
    ranges.emplace_back(0, 0, false, 0, 0, false, 0, 0);
    ranges.emplace_back(0, fw_seed_1.size(), false,
                        0, fw_seed_1.size(), false,
                        0, 0);
    ranges.emplace_back(fw_seed_1.size() + target_1_truth.size(), fw_seed_1.size() + target_1_truth.size() + rc_seed_target.size(), false,
                        fw_seed_1.size(), fw_seed_1.size() + rc_seed_query.size(), true,
                        0, 0);
    ranges.emplace_back(fw_seed_1.size() + target_1_truth.size() + rc_seed_target.size() + target_2_truth.size(), target.size(), false,
                        fw_seed_1.size() + rc_seed_query.size() + query_1_truth.size() + query_rc_2_truth.size(), query.size(), false,
                        0, 0);
    ranges.emplace_back(target.size(), target.size(), false, query.size(), query.size(), false, 0, 0);

    ScoreModel score_model;

    auto [ql_1, query_w_1, query_rc_w_1, target_w_1, mum_length_1, qi_1, ti_1, eml_1]
        = extract_gap_seqs_switch<false>(target, query, query_rc, ranges, 1, 2, 2, 3);
    ASSERT_FALSE(ql_1);

    auto [ql_2, query_w_2, query_rc_w_2, target_w_2, mum_length_2, qi_2, ti_2, eml_2]
        = extract_gap_seqs_switch<true>(target, query, query_rc, ranges, 1, 2, 2, 3);
    ASSERT_TRUE(ql_2);

    ASSERT_EQ(query_1_truth, query_w_1);
    ASSERT_EQ(query_rc_1_truth, query_rc_w_1);
    ASSERT_EQ(target_1_truth, target_w_1);
    ASSERT_EQ(query_2_truth, query_w_2);
    ASSERT_EQ(query_rc_2_truth, query_rc_w_2);
    ASSERT_EQ(target_2_truth, target_w_2);

    auto aligner = make_aligner(score_model);

    auto [score_1, cigar_1, r_consumed_1, q_consumed_1, inv_length_1, inv_length_r_1,
          score_2, cigar_2, r_consumed_2, q_consumed_2, inv_length_2, inv_length_r_2,
          inv_len, inv_len_r]
        = run_alignment(*aligner, score_model,
                        query, query_rc, target,
                        query_w_1, query_rc_w_1, target_w_1,
                        query_w_2, query_rc_w_2, target_w_2);

    EXPECT_EQ(Cigar("3D15i"), cigar_1);
    EXPECT_EQ(Cigar("8="), cigar_2);

    EXPECT_EQ(15, inv_len);
    EXPECT_EQ(15, inv_len_r);
}

TEST(CIGARTest, NoN) {
    std::string target = "ATGC";
    std::string query = "ATGG";
    Cigar cigar("3=1X");
    EXPECT_EQ(cigar, cigar_fix_n(cigar, target, query));
}