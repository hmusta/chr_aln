#include <gtest/gtest.h>

#include "helpers.hpp"

#include <string>
#include <string_view>

#include "chaining.hpp"

TEST(ChainingTest, NoOverlapNoInv) {
    std::string shared_prefix = "TTTAGGATTGCAGATGATTTTGCTTTTAGTTATATTTTCCACACACTGTCCAGAATTACC";
    std::string target_middle = "AGATGTTATACACAAAATGTCGAGCAGACTTCACTTGAGAAAACATTCTTTCTCAGGATT";
    std::string shared_suffix = "CCAGTCACAACCTGCAATTGTGCGGCAAAGTGCAGCAAAGTATTTCCTAAATCTGCAAAC";
    std::string query_middle = reverse_complement(target_middle);

    std::string target = shared_prefix + target_middle + shared_suffix;
    std::string query = shared_prefix + query_middle + shared_suffix;
    std::string query_rc = reverse_complement(query);

    ASSERT_EQ(target.size(), query.size());

    std::vector<Ranges> ranges;
    ranges.emplace_back(0, shared_prefix.size(), false,
                        0, shared_prefix.size(), false,
                        0, 0);
    ranges.emplace_back(shared_prefix.size(), shared_prefix.size() + target_middle.size(), false,
                        shared_prefix.size(), shared_prefix.size() + target_middle.size(), true,
                        0, 0);
    ranges.emplace_back(shared_prefix.size() + target_middle.size(), target.size(), false,
                        shared_prefix.size() + target_middle.size(), query.size(), false,
                        0, 0);

    ScoreModel score_model(1, -9, -16, -2, -41, -1, -41, 0);
    double exp_mismatch_frac_between_mum_bp = 0.01;
    SOffset max_gap = max_offset;
    ChainScoreModel chain_score_model(target_middle.size(), score_model.gap_switch, exp_mismatch_frac_between_mum_bp, max_gap);

    auto [chain, chain_score]
            = chain_ranges(target, query, query_rc, ranges, score_model, chain_score_model);

    ASSERT_LE(2u, chain.size());
    EXPECT_EQ(0u, chain.front().size());
    EXPECT_EQ(0u, chain.back().size());

    EXPECT_EQ(ranges[0], chain[1]);

    ASSERT_LT(2u, chain.size());
    EXPECT_EQ(ranges[2], chain[2]);

    EXPECT_EQ(4u, chain.size());
}

TEST(ChainingTest, NoOverlap) {
    std::string shared_prefix = "TTTAGGATTGCAGATGATTTTGCTTTTAGTTATATTTTCCACACACTGTCCAGAATTACC";
    std::string target_middle = "AGATGTTATACACAAAATGTCGAGCAGACTTCACTTGAGAAAACATTCTTTCTCAGGATT";
    std::string shared_suffix = "CCAGTCACAACCTGCAATTGTGCGGCAAAGTGCAGCAAAGTATTTCCTAAATCTGCAAAC";
    std::string query_middle = reverse_complement(target_middle);

    std::string target = shared_prefix + target_middle + shared_suffix;
    std::string query = shared_prefix + query_middle + shared_suffix;
    std::string query_rc = reverse_complement(query);

    ASSERT_EQ(target.size(), query.size());

    std::vector<Ranges> ranges;
    ranges.emplace_back(0, shared_prefix.size(), false,
                        0, shared_prefix.size(), false,
                        0, 0);
    ranges.emplace_back(shared_prefix.size(), shared_prefix.size() + target_middle.size(), false,
                        shared_prefix.size(), shared_prefix.size() + target_middle.size(), true,
                        0, 0);
    ranges.emplace_back(shared_prefix.size() + target_middle.size(), target.size(), false,
                        shared_prefix.size() + target_middle.size(), query.size(), false,
                        0, 0);

    ScoreModel score_model(1, -9, -16, -2, -41, -1, -41, 0);
    double exp_mismatch_frac_between_mum_bp = 0.01;
    SOffset max_gap = max_offset;
    ChainScoreModel chain_score_model(target_middle.size(), score_model.gap_switch, exp_mismatch_frac_between_mum_bp, max_gap);

    auto [chain, chain_score]
            = chain_ranges(target, query, query_rc, ranges, score_model, chain_score_model, true);

    ASSERT_LE(2u, chain.size());
    EXPECT_EQ(0u, chain.front().size());
    EXPECT_EQ(0u, chain.back().size());

    EXPECT_EQ(ranges[0], chain[1]);

    ASSERT_LT(2u, chain.size());
    EXPECT_EQ(ranges[1], chain[2]);

    ASSERT_LT(3u, chain.size());
    EXPECT_EQ(ranges[2], chain[3]);

    EXPECT_EQ(5u, chain.size());
}

TEST(ChainingTest, Overlap) {
    std::string shared_prefix = "TTTAGGATTGCAGATGATTTTGCTTTTAGTTATATTTTCCACACACTGTCCAGAATTACC";
    std::string target_middle = "AGATGTTATACACAAAATGTCGAGCAGACTTCACTTGAGAAAACATTCTTTCTCAGGATT";
    std::string shared_suffix = "CCAGTCACAACCTGCAATTGTGCGGCAAAGTGCAGCAAAGTATTTCCTAAATCTGCAAAC";
    std::string query_middle = reverse_complement(target_middle);

    std::string target = shared_prefix + target_middle + shared_suffix;
    std::string query = shared_prefix + query_middle + shared_suffix;
    std::string query_rc = reverse_complement(query);

    ASSERT_EQ(target.size(), query.size());

    std::vector<Ranges> ranges;
    ranges.emplace_back(0, shared_prefix.size() + 1, false,
                        0, shared_prefix.size() + 1, false,
                        0, 0);
    ranges.emplace_back(shared_prefix.size(), shared_prefix.size() + target_middle.size() + 1, false,
                        shared_prefix.size(), shared_prefix.size() + target_middle.size() + 1, true,
                        0, 0);
    ranges.emplace_back(shared_prefix.size() + target_middle.size() - 1, target.size(), false,
                        shared_prefix.size() + target_middle.size() - 1, query.size(), false,
                        0, 0);

    ScoreModel score_model(1, -9, -16, -2, -41, -1, -41, 0);
    double exp_mismatch_frac_between_mum_bp = 0.01;
    SOffset max_gap = max_offset;
    ChainScoreModel chain_score_model(target_middle.size(), score_model.gap_switch, exp_mismatch_frac_between_mum_bp, max_gap);

    auto [chain, chain_score]
            = chain_ranges(target, query, query_rc, ranges, score_model, chain_score_model, true);

    ASSERT_LE(2u, chain.size());
    EXPECT_EQ(0u, chain.front().size());
    EXPECT_EQ(0u, chain.back().size());

    EXPECT_EQ(ranges[0], chain[1]);

    ASSERT_LT(2u, chain.size());
    EXPECT_EQ(ranges[1], chain[2]);

    ASSERT_LT(3u, chain.size());
    EXPECT_EQ(ranges[2], chain[3]);

    EXPECT_EQ(5u, chain.size());
}