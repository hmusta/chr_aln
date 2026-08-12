#include <gtest/gtest.h>

#include "helpers.hpp"

#include <sstream>
#include <string>
#include <vector>

#include "monotone_diagonal_profile.hpp"
#include "offset_vector.hpp"
#include "readers.hpp"

// ------------------------------------------------------------ OffsetVector --

TEST(OffsetVectorTest, DeallocFrontKeepsLogicalIndices) {
    OffsetVector<int> v;
    v.emplace_back(10);
    v.emplace_back(11);
    v.emplace_back(12);

    EXPECT_EQ(3u, v.size());
    EXPECT_EQ(3u, v.alloc_size());
    EXPECT_EQ(0u, v.offset());
    EXPECT_EQ(10, v[0]);

    v.dealloc_front(2);

    // indices stay logical: v[2] is still the third element pushed
    EXPECT_EQ(3u, v.size());
    EXPECT_EQ(1u, v.alloc_size());
    EXPECT_EQ(2u, v.offset());
    EXPECT_EQ(12, v[2]);
}

// clear() used to drop the elements without resetting the offset, so an
// emptied vector still reported a non-zero size and mis-indexed everything.
TEST(OffsetVectorTest, ClearResetsOffset) {
    OffsetVector<int> v;
    v.emplace_back(1);
    v.emplace_back(2);
    v.emplace_back(3);
    v.dealloc_front(2);
    ASSERT_EQ(2u, v.offset());

    v.clear();

    EXPECT_EQ(0u, v.size());
    EXPECT_EQ(0u, v.alloc_size());
    EXPECT_EQ(0u, v.offset());
}

TEST(OffsetVectorTest, ResizeGrowsPastTheOffset) {
    OffsetVector<int> v;
    v.emplace_back(1);
    v.emplace_back(2);
    v.dealloc_front();
    ASSERT_EQ(1u, v.offset());

    v.resize(5);
    EXPECT_EQ(5u, v.size());
    EXPECT_EQ(4u, v.alloc_size());
    EXPECT_EQ(2, v[1]);
}

// -------------------------------------------------------------- ScoreModel --

TEST(ScoreModelTest, DefaultIsEditDistance) {
    ScoreModel model;
    EXPECT_EQ("Edit distance", model.get_model_str());
}

TEST(ScoreModelTest, ZeroGapOpenIsLinear) {
    // o1 == o2 == 0 and e1 == e2: no affine component
    ScoreModel model(2, -4, 0, -2, 0, -2, 0, 0);
    EXPECT_EQ("Gap linear", model.get_model_str());
}

// A single-piece affine model used to be reported as "Gap linear", because the
// classifier tested `gap_open_s > 0` while the constructor asserts it is <= 0.
TEST(ScoreModelTest, NegativeGapOpenIsAffine) {
    ScoreModel model(2, -4, -4, -2, -4, -2, 0, 0);
    ASSERT_LT(model.gap_open_s, 0);
    EXPECT_EQ("Gap affine", model.get_model_str());
}

TEST(ScoreModelTest, DifferingPiecesAreTwoPieceAffine) {
    ScoreModel model(2, -4, -4, -2, -8, -1, 0, 0);
    EXPECT_EQ("2-piece gap affine", model.get_model_str());

    // the long-gap piece must be cheaper beyond the switch length
    ASSERT_GT(model.gap_switch, 1);
    ASSERT_LT(model.gap_switch, max_offset);
    EXPECT_LE(model.get_gap_penalty(model.gap_switch),
              model.gap_open_p + model.gap_ext_p * model.gap_switch);
}

TEST(ScoreModelTest, GapPenaltyTakesTheCheaperPiece) {
    ScoreModel model(2, -4, -4, -2, -8, -1, 0, 0);
    for (SOffset len = 1; len < 20; ++len) {
        EXPECT_EQ(std::min(model.gap_open_p + model.gap_ext_p * len,
                           model.gap_open2_p + model.gap_ext2_p * len),
                  model.get_gap_penalty(len)) << "at length " << len;
    }
    EXPECT_EQ(0u, model.get_gap_penalty(0));
}

// ---------------------------------------------------------------- sequence --

TEST(SequenceTest, SanitizeNucUpperCases) {
    EXPECT_EQ('A', sanitize_nuc('a'));
    EXPECT_EQ('C', sanitize_nuc('c'));
    EXPECT_EQ('G', sanitize_nuc('g'));
    EXPECT_EQ('T', sanitize_nuc('t'));
    EXPECT_EQ('N', sanitize_nuc('n'));
}

// Bytes with the high bit set used to index seq_ntext_table at a negative
// offset, because the char was cast to a signed type.
TEST(SequenceTest, SanitizeNucMapsHighBytesToN) {
    for (int byte : { 0x80, 0xA5, 0xC3, 0xFF }) {
        EXPECT_EQ('N', sanitize_nuc(static_cast<char>(byte)))
            << "byte 0x" << std::hex << byte;
    }
}

TEST(SequenceTest, ReverseComplementRoundTrips) {
    const std::string fw = "ACGTACGGTTCA";
    const std::string rc = reverse_complement(fw);

    EXPECT_EQ(fw.size(), rc.size());
    EXPECT_EQ("TGAACCGTACGT", rc);
    EXPECT_EQ(fw, reverse_complement(rc));
    EXPECT_TRUE(is_reverse_complement(fw, rc));
}

TEST(SequenceTest, ReverseComplementHandlesAmbiguityCodes) {
    // N is its own complement, R<->Y, K<->M, B<->V, D<->H, S and W are self
    EXPECT_EQ("N", reverse_complement("N"));
    EXPECT_EQ("Y", reverse_complement("R"));
    EXPECT_EQ("R", reverse_complement("Y"));
    EXPECT_EQ("M", reverse_complement("K"));
    EXPECT_EQ("S", reverse_complement("S"));
    EXPECT_EQ("W", reverse_complement("W"));
}

// The complement table is indexed by the character, so a byte >= 0x80 used to
// read out of bounds. Such bytes have no complement and map to themselves.
TEST(SequenceTest, ReverseComplementHandlesHighBytes) {
    std::string fw = "ACGT";
    fw.push_back(static_cast<char>(0xC3));

    const std::string rc = reverse_complement(fw);
    ASSERT_EQ(fw.size(), rc.size());
    EXPECT_EQ(static_cast<char>(0xC3), rc.front());
    EXPECT_EQ("ACGT", rc.substr(1));
}

TEST(SequenceTest, IsReverseComplementRejectsMismatches) {
    // ACGT is a palindrome, so it really is its own reverse complement
    EXPECT_TRUE(is_reverse_complement("ACGT", "ACGT"));

    // ACGG is not
    EXPECT_FALSE(is_reverse_complement("ACGG", "ACGG"));
    EXPECT_TRUE(is_reverse_complement("ACGG", "CCGT"));

    // differing lengths never match
    EXPECT_FALSE(is_reverse_complement("ACGT", "ACG"));
}

TEST(SequenceTest, HasLargeGapCountsNsPerSequence) {
    const std::string clean(100, 'A');
    const std::string gappy(100, 'N');

    EXPECT_FALSE(has_large_gap(clean, clean, 10));
    EXPECT_TRUE(has_large_gap(gappy, clean, 10));
    EXPECT_TRUE(has_large_gap(clean, gappy, 10));
}

// ------------------------------------------------------------------ reader --

TEST(ReaderTest, ReadFastaUpperCasesAndSanitizes) {
    std::istringstream fin(">chr1 some description\nacgt\nACGT\n");
    auto [seq, header] = read_fasta(fin);

    EXPECT_EQ("chr1", header);
    EXPECT_EQ("ACGTACGT", seq);
}

TEST(ReaderTest, ReadFastaReplacesUnknownBytesWithN) {
    std::string fasta = ">chr1\nAC";
    fasta.push_back(static_cast<char>(0xC3));
    fasta += "gt\n";

    std::istringstream fin(fasta);
    auto [seq, header] = read_fasta(fin);

    EXPECT_EQ("chr1", header);
    EXPECT_EQ("ACNGT", seq);
}

// --------------------------------------------------- MonotoneDiagonalProfile --

TEST(DiagonalProfileTest, ForwardStoresDisjointAscendingRecords) {
    chr_aln::ForwardDiagonalProfile profile;
    EXPECT_TRUE(profile.empty());
    EXPECT_FALSE(profile.get_global_min().has_value());

    EXPECT_TRUE(profile.extend(0, 10, Penalty(5)));
    EXPECT_TRUE(profile.extend(10, 20, Penalty(7)));

    ASSERT_EQ(2u, profile.size());
    EXPECT_EQ(Penalty(5), profile.get_global_min().value());
    EXPECT_EQ(0, profile.get_global_begin());
    EXPECT_EQ(20, profile.get_global_end());
}

TEST(DiagonalProfileTest, ForwardExtendsAnEqualPenaltyRun) {
    chr_aln::ForwardDiagonalProfile profile;
    ASSERT_TRUE(profile.extend(0, 10, Penalty(5)));

    // same penalty and a longer reach: the existing record grows
    EXPECT_TRUE(profile.extend(0, 15, Penalty(5)));
    EXPECT_EQ(1u, profile.size());
    EXPECT_EQ(15, profile.get_global_end());

    // an empty or non-advancing range changes nothing
    EXPECT_FALSE(profile.extend(5, 5, Penalty(5)));
    EXPECT_FALSE(profile.extend(0, 12, Penalty(5)));
    EXPECT_EQ(1u, profile.size());
}

TEST(DiagonalProfileTest, ForwardQueryFindsTheCoveringRecord) {
    chr_aln::ForwardDiagonalProfile profile;
    ASSERT_TRUE(profile.extend(0, 10, Penalty(5)));
    ASSERT_TRUE(profile.extend(10, 20, Penalty(7)));

    auto [begin, end] = profile.query_min(0, 5);
    ASSERT_NE(begin, end);
    EXPECT_EQ(Penalty(5), begin->penalty);

    auto [begin2, end2] = profile.query_min(12, 15);
    ASSERT_NE(begin2, end2);
    EXPECT_EQ(Penalty(7), begin2->penalty);

    // a query past everything stored comes back empty
    auto [begin3, end3] = profile.query_min(100, 200);
    EXPECT_EQ(begin3, end3);
}

// The descending specialisation is never instantiated by the aligner, so a
// compile error in its branch went unnoticed. Instantiating it here keeps that
// branch built.
TEST(DiagonalProfileTest, BackwardProfileInstantiates) {
    chr_aln::BackwardDiagonalProfile profile;
    EXPECT_TRUE(profile.empty());

    auto [begin, end] = profile.query_min(0, 10);
    EXPECT_EQ(begin, end);

    auto [rbegin, rend] = profile.query_min_range(0, 10);
    EXPECT_EQ(rbegin, rend);
}
// --------------------------------------------------------------- Ranges --
TEST(RangesTest, TrimAccumulatesOnBothEnds) {
    Ranges r(100, 200, false, 300, 400, false);
    EXPECT_EQ(0, r.left_trim);
    EXPECT_EQ(0, r.right_trim);

    r.trim_prefix(7);
    r.trim_suffix(3);

    EXPECT_EQ(7, r.left_trim);
    EXPECT_EQ(3, r.right_trim);
    EXPECT_EQ(107, r.rbegin);
    EXPECT_EQ(197, r.rend);
    EXPECT_EQ(307, r.qbegin);
    EXPECT_EQ(397, r.qend);
    EXPECT_EQ(90u, r.size());
}

TEST(RangesTest, TrimOnAnInvertedRangeMovesTheOppositeQueryEnd) {
    Ranges r(100, 200, false, 300, 400, true);

    r.trim_prefix(10);
    EXPECT_EQ(110, r.rbegin);
    EXPECT_EQ(390, r.qend);     // inverted: the query shrinks from the far end
    EXPECT_EQ(90u, r.size());

    r.trim_suffix(5);
    EXPECT_EQ(195, r.rend);
    EXPECT_EQ(305, r.qbegin);
    EXPECT_EQ(85u, r.size());
}

// The constructor's trailing parameters were unnamed, so the trim counts they
// were given were silently replaced with zero.
TEST(RangesTest, ConstructorKeepsTheTrimCounts) {
    Ranges r(100, 200, false, 300, 400, false, 7, 3);
    EXPECT_EQ(7, r.left_trim);
    EXPECT_EQ(3, r.right_trim);
}

TEST(RangesTest, CopyKeepsTheTrimCounts) {
    Ranges r(100, 200, false, 300, 400, false);
    r.trim_prefix(7);
    r.trim_suffix(3);

    const Ranges copy = r;
    EXPECT_EQ(7, copy.left_trim);
    EXPECT_EQ(3, copy.right_trim);
}

TEST(RangesTest, RcRangesFlipsTheQueryAndKeepsTheTrims) {
    const SOffset query_size = 1000;

    Ranges r(100, 200, false, 300, 400, false);
    r.trim_prefix(7);
    r.trim_suffix(3);

    const std::vector<Ranges> flipped = rc_ranges({ r }, query_size);
    ASSERT_EQ(1u, flipped.size());

    // reference coordinates are untouched, query coordinates mirror
    EXPECT_EQ(r.rbegin, flipped[0].rbegin);
    EXPECT_EQ(r.rend, flipped[0].rend);
    EXPECT_EQ(query_size - r.qend, flipped[0].qbegin);
    EXPECT_EQ(query_size - r.qbegin, flipped[0].qend);
    EXPECT_NE(r.qorientation, flipped[0].qorientation);

    EXPECT_EQ(r.left_trim, flipped[0].left_trim);
    EXPECT_EQ(r.right_trim, flipped[0].right_trim);
}

TEST(RangesTest, RcRangesRoundTrips) {
    const SOffset query_size = 1000;
    const Ranges r(100, 200, false, 300, 400, false, 7, 3);

    const std::vector<Ranges> once = rc_ranges({ r }, query_size);
    const std::vector<Ranges> twice = rc_ranges(once, query_size);

    ASSERT_EQ(1u, twice.size());
    EXPECT_EQ(r, twice[0]);
    EXPECT_EQ(r.left_trim, twice[0].left_trim);
    EXPECT_EQ(r.right_trim, twice[0].right_trim);
}

// ------------------------------------------------------------ read_mummer --

namespace {

// target and query share "ACGTACGT" at target[0..8) and query[2..10)
struct MummerFixture {
    std::string target = "ACGTACGTAA";
    std::string query = "TTACGTACGT";
    std::string query_rc = reverse_complement(query);
    std::string theader = "ref";
    std::string qheader = "qry";

    std::vector<Ranges> read(const std::string &text, SOffset min_len = 8) {
        std::istringstream fin(text);
        return read_mummer(fin, target, theader, query, qheader, query_rc, min_len);
    }
};

} // namespace

TEST(ReadMummerTest, ParsesAForwardMatch) {
    MummerFixture f;
    const std::vector<Ranges> mums = f.read("> qry\n1 3 8\n");

    ASSERT_EQ(1u, mums.size());
    EXPECT_EQ(0, mums[0].rbegin);
    EXPECT_EQ(8, mums[0].rend);
    EXPECT_EQ(2, mums[0].qbegin);
    EXPECT_EQ(10, mums[0].qend);
    EXPECT_FALSE(mums[0].qorientation);
}

// MUMmer output ends with a newline, so the final getline hands back an empty
// line. It used to be parsed as a record, leaving the coordinates untouched.
TEST(ReadMummerTest, SkipsBlankTrailingLines) {
    MummerFixture f;
    EXPECT_EQ(1u, f.read("> qry\n1 3 8\n").size());
    EXPECT_EQ(1u, f.read("> qry\n1 3 8\n\n").size());
    EXPECT_EQ(1u, f.read("> qry\n1 3 8\n   \n\t\n").size());
}

TEST(ReadMummerTest, RejectsATruncatedLine) {
    MummerFixture f;
    EXPECT_THROW(f.read("> qry\n1 3\n"), std::runtime_error);
    EXPECT_THROW(f.read("> qry\nnot a mum\n"), std::runtime_error);
}

// A negative field parses "successfully" into an unsigned type and wraps to a
// huge value, so the coordinates have to be read signed and range-checked.
TEST(ReadMummerTest, RejectsNegativeCoordinates) {
    MummerFixture f;
    EXPECT_THROW(f.read("> qry\n-1 3 8\n"), std::runtime_error);
    EXPECT_THROW(f.read("> qry\n1 -3 8\n"), std::runtime_error);
    EXPECT_THROW(f.read("> qry\n1 3 -8\n"), std::runtime_error);
    EXPECT_THROW(f.read("> qry\n0 3 8\n"), std::runtime_error);
}

TEST(ReadMummerTest, IgnoresExtraColumns) {
    MummerFixture f;
    // some MUMmer modes emit a fourth field
    EXPECT_EQ(1u, f.read("> qry\n1 3 8 extra\n").size());
}
