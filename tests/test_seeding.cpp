#include <gtest/gtest.h>

#include "repeat_aligner.hpp"

#include <string>
#include <string_view>

TEST(SeedTest, CallMAMs) {
    std::string target = "ATGC";
    std::string query = "ATGC";
    std::string query_rc = "GCAT";

    std::vector<Ranges> mums;
    call_mums(query, query_rc, target, 4, [&mums](Ranges&& mum) {
        mums.emplace_back(std::move(mum));
    });

    ASSERT_EQ(1u, mums.size());
    EXPECT_EQ(mums[0], Ranges(0, 4, false, 0, 4, false));
}

TEST(SeedTest, CallMAMsWithN) {
    std::string target = "ATGCN";
    std::string query = "ATGCN";
    std::string query_rc = "NGCAT";

    std::vector<Ranges> mums;
    call_mums(query, query_rc, target, 4, [&mums](Ranges&& mum) {
        mums.emplace_back(std::move(mum));
    });

    ASSERT_EQ(1u, mums.size());
    EXPECT_EQ(mums[0], Ranges(0, 4, false, 0, 4, false));
}