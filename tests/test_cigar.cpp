#include <gtest/gtest.h>

#include "helpers.hpp"

#include <sstream>
#include <string>
#include <vector>

#include "cigar.hpp"

namespace {

Cigar make_cigar(const std::vector<std::pair<char, Offset>> &ops) {
    Cigar cigar;
    for (const auto &[op, num] : ops)
        cigar.push(op, num);
    return cigar;
}

std::string to_string(const Cigar &cigar) {
    std::ostringstream sout;
    sout << cigar;
    return sout.str();
}

// total length of every operation, which merging must never change
Offset total_length(const Cigar &cigar) {
    Offset total = 0;
    for (const auto &[op, num] : cigar)
        total += num;
    return total;
}

} // namespace

TEST(CigarTest, InsertWithoutMerge) {
    Cigar dst = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src = make_cigar({ { EQ_OP, 7 } });

    dst.insert(dst.begin() + 1, src.begin(), src.end());

    EXPECT_EQ("10=7=5X", to_string(dst));
    EXPECT_EQ(22u, total_length(dst));
}

TEST(CigarTest, InsertMergesWithFollowingOp) {
    Cigar dst = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src = make_cigar({ { NEQ_OP, 7 } });

    // the range's last op matches the op at the insertion point
    dst.insert(dst.begin() + 1, src.begin(), src.end(), true);

    EXPECT_EQ("10=12X", to_string(dst));
    EXPECT_EQ(22u, total_length(dst));
}

TEST(CigarTest, InsertMergesWithPrecedingOp) {
    Cigar dst = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src = make_cigar({ { EQ_OP, 7 } });

    // the range's first op matches the op before the insertion point
    dst.insert(dst.begin() + 1, src.begin(), src.end(), true);

    EXPECT_EQ("17=5X", to_string(dst));
    EXPECT_EQ(22u, total_length(dst));
}

TEST(CigarTest, InsertMergesOnBothSides) {
    Cigar dst = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src = make_cigar({ { EQ_OP, 7 }, { NEQ_OP, 3 } });

    dst.insert(dst.begin() + 1, src.begin(), src.end(), true);

    EXPECT_EQ("17=8X", to_string(dst));
    EXPECT_EQ(25u, total_length(dst));
}

// A one-element range whose op matches both neighbours. Merging it into the
// following op empties the range, and it must not then also be merged into the
// preceding one: that double-counted the run and left begin past end, which
// std::vector::insert rejects with a length_error.
TEST(CigarTest, InsertSingleOpBetweenTwoMatchingOps) {
    Cigar dst = make_cigar({ { NEQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src = make_cigar({ { NEQ_OP, 7 } });

    dst.insert(dst.begin() + 1, src.begin(), src.end(), true);

    EXPECT_EQ(22u, total_length(dst));
    EXPECT_EQ("10X12X", to_string(dst));
}

TEST(CigarTest, InsertSingleOpMatchingPrecedingOnly) {
    Cigar dst = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src = make_cigar({ { EQ_OP, 7 } });

    // insertion point is the end, so only the preceding op can merge
    dst.insert(dst.end(), src.begin(), src.end(), true);

    EXPECT_EQ("10=5X7=", to_string(dst));
    EXPECT_EQ(22u, total_length(dst));
}

TEST(CigarTest, InsertAtBeginMergesWithFollowingOnly) {
    Cigar dst = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src = make_cigar({ { EQ_OP, 7 } });

    dst.insert(dst.begin(), src.begin(), src.end(), true);

    EXPECT_EQ("17=5X", to_string(dst));
    EXPECT_EQ(22u, total_length(dst));
}

TEST(CigarTest, InsertEmptyRangeIsANoOp) {
    Cigar dst = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar src;

    dst.insert(dst.begin() + 1, src.begin(), src.end(), true);

    EXPECT_EQ("10=5X", to_string(dst));
    EXPECT_EQ(15u, total_length(dst));
}

TEST(CigarTest, InsertIntoEmptyCigar) {
    Cigar dst;
    Cigar src = make_cigar({ { EQ_OP, 7 }, { NEQ_OP, 3 } });

    dst.insert(dst.begin(), src.begin(), src.end(), true);

    EXPECT_EQ("7=3X", to_string(dst));
    EXPECT_EQ(10u, total_length(dst));
}

// Merging must leave the alignment itself untouched: same operations, same
// consumed lengths, only the run boundaries differ.
TEST(CigarTest, MergePreservesConsumedLengths) {
    Cigar merged = make_cigar({ { EQ_OP, 10 }, { NEQ_OP, 5 } });
    Cigar plain = merged;
    Cigar src = make_cigar({ { NEQ_OP, 7 }, { TARGET_CONSUME_OP, 4 } });

    merged.insert(merged.begin() + 1, src.begin(), src.end(), true);
    plain.insert(plain.begin() + 1, src.begin(), src.end());

    auto consumed = [](const Cigar &cigar) {
        Offset target = 0;
        Offset query = 0;
        for (const auto &[op, num] : cigar) {
            target += consumes_target(op) ? num : 0;
            query += consumes_query(op) ? num : 0;
        }
        return std::make_pair(target, query);
    };

    EXPECT_EQ(consumed(plain), consumed(merged));
    EXPECT_EQ(total_length(plain), total_length(merged));
    EXPECT_LE(merged.data_.size(), plain.data_.size());
}
// cigar_fix_n relabels aligned positions against the sequences: = and X only
// ever cover positions where neither side has an N, and M covers the rest.
TEST(CigarTest, FixNSplitsMatchesAroundNs) {
    const std::string target = "ACGTNNAC";
    const std::string query  = "ACGTNNAC";

    EXPECT_EQ("4=2M2=", to_string(cigar_fix_n(make_cigar({ { EQ_OP, 8 } }), target, query)));
}

TEST(CigarTest, FixNSeparatesMatchesFromMismatches) {
    const std::string target = "ACGTACGT";
    const std::string query  = "ACGAACGT";

    // one substitution at offset 3
    EXPECT_EQ("3=1X4=", to_string(cigar_fix_n(make_cigar({ { EQ_OP, 8 } }), target, query)));
}

TEST(CigarTest, FixNTreatsAnNOnEitherSideAsM) {
    const std::string target = "ACNNAC";
    const std::string query  = "ACGTAC";   // Ns only on the target side

    EXPECT_EQ("2=2M2=", to_string(cigar_fix_n(make_cigar({ { EQ_OP, 2 }, { NEQ_OP, 2 }, { EQ_OP, 2 } }),
                                              target, query)));
}

TEST(CigarTest, FixNLeavesIndelsAlone) {
    const std::string target = "ACGTAC";
    const std::string query  = "ACGAC";

    const Cigar fixed = cigar_fix_n(
        make_cigar({ { EQ_OP, 3 }, { TARGET_CONSUME_OP, 1 }, { EQ_OP, 2 } }), target, query);
    EXPECT_EQ("3=1D2=", to_string(fixed));
}

TEST(CigarTest, FixNIsIdempotent) {
    const std::string target = "ACNNACGTAC";
    const std::string query  = "ACGTACGAAC";

    const Cigar once = cigar_fix_n(make_cigar({ { EQ_OP, 10 } }), target, query);
    const Cigar twice = cigar_fix_n(once, target, query);
    EXPECT_EQ(to_string(once), to_string(twice));
}
