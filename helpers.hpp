#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>


using Offset = uint64_t;
using SOffset = std::make_signed_t<Offset>;

static constexpr SOffset max_offset = std::numeric_limits<SOffset>::max() / 4;

// this needs to be -1 for the max operations in wfa_extend and wfa_next to work
static constexpr Offset npos = -1;

using Penalty = uint64_t;
using Score = std::make_signed_t<Penalty>;

// ensure that the asserts for positiveness below work
static_assert(static_cast<Score>(static_cast<Penalty>(SOffset(-1))) == SOffset(-1));

using Diag = int32_t;
constexpr Diag min_diag = std::numeric_limits<Diag>::min() / 4;
constexpr Diag max_diag = std::numeric_limits<Diag>::max() / 4;
constexpr Diag max_diag_width = max_diag - min_diag;

using UDiag = std::make_unsigned_t<Diag>;

using SeqPair = std::pair<std::string_view, std::string_view>;

static const unsigned char alph_size_ext = 15;
static const unsigned char seq_ntext_table[256]
    = { 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  0, 10,  1, 12, 15, 15,  2,
        13, 15, 15,  8, 15,  9, 14, 15, 15, 15,  6,  4,  3,  3, 11,  5, 15,  7, 15, 15, 15, 15, 15, 15,
        15,  0, 10,  1, 12, 15, 15,  2, 13, 15, 15,  8, 15,  9, 14, 15, 15, 15,  6,  4,  3,  3, 11,  5,
        15,  7, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };

inline char sanitize_nuc(unsigned char c) {
    c = std::toupper(c);
    assert(seq_ntext_table[c] <= alph_size_ext);
    if (seq_ntext_table[c] >= alph_size_ext)
        c = 'N';

    assert(seq_ntext_table[c] < alph_size_ext);
    return c;
}

std::string reverse_complement(std::string_view fw);
bool is_reverse_complement(std::string_view fw, std::string_view rc);

inline char orientation_to_char(bool orientation) {
    static const char tab[] = { '+', '-' };
    return tab[orientation];
}

struct nuc_cmp_n_eq {
    bool operator()(unsigned char a, unsigned char b) const {
        assert(seq_ntext_table[a] < alph_size_ext);
        assert(seq_ntext_table[b] < alph_size_ext);
        return a == b || a == 'N' || b == 'N';
    }
};

template <typename T, typename U>
inline T get_mismatch_it(T q_begin, T q_end, U t_begin, U t_end) {
    auto it = q_begin;
    auto jt = t_begin;
    while (it != q_end && jt != t_end) {
        std::tie(it, jt) = std::mismatch(it, q_end, jt, t_end);
        bool found = false;
        while (it != q_end && jt != t_end && (*it == 'N' || *jt == 'N')) {
            found = true;
            ++it;
            ++jt;
        }

        if (!found)
            break;
    }

    return it;
}

inline int match_char(int v, int h, void* stringview_pair) {
    const auto& [q, t] = *reinterpret_cast<const SeqPair*>(stringview_pair);

    int qlen = q.size();
    int rlen = t.size();

    // Check boundaries
    if (v >= qlen || h >= rlen)
        return 0;

    assert(v < qlen && h < rlen);

    // Compare arrays
    return q[v] == t[h] || q[v] == 'N' || t[h] == 'N';
}

inline bool has_large_gap(std::string_view query,
                          std::string_view target,
                          size_t cutoff = 10000) {
    size_t n_count = 0;
    for (char c : query) {
        if (c == 'N') {
            ++n_count;
            if (n_count >= cutoff)
                break;
        }
    }

    if (n_count >= cutoff)
        return true;

    n_count = 0;
    for (char c : target) {
        if (c == 'N') {
            ++n_count;
            if (n_count >= cutoff)
                break;
        }
    }

    return n_count >= cutoff;
}

struct Ranges {
    Ranges(SOffset rb, SOffset re, bool ro, SOffset qb, SOffset qe, bool qo, SOffset lt = 0, SOffset rt = 0)
        : rbegin(rb), rend(re), rorientation(ro),
          qbegin(qb), qend(qe), qorientation(qo),
          left_trim(lt), right_trim(rt) { assert_valid(); }

    static_assert(std::is_same_v<SOffset, std::make_signed_t<SOffset>>);
    SOffset rbegin;
    SOffset rend;
    bool rorientation;
    SOffset qbegin;
    SOffset qend;
    bool qorientation;
    SOffset left_trim;
    SOffset right_trim;

    size_t size() const {
        assert_valid();
        return rend - rbegin;
    }

    Diag diag() const {
        assert_valid();
        return rend - qend;
    }

    void assert_valid(bool allow_empty = true) const {
        assert(!rorientation);
        assert(rbegin >= 0);
        assert(qbegin >= 0);
        assert(rbegin <= rend);
        assert(allow_empty || rbegin < rend);
        assert(qbegin <= qend);
        assert(allow_empty || qbegin < qend);
        assert(rend - rbegin == qend - qbegin);

        #ifdef NDEBUG
        std::ignore = allow_empty;
        #endif
    }

    void shift_start(size_t target_shift, size_t query_shift) {
        assert_valid();
        rbegin += target_shift;
        rend += target_shift;
        qbegin += query_shift;
        qend += query_shift;
        assert_valid();
    }

    void trim_prefix(size_t trim, bool allow_empty = true) {
        assert_valid();
        assert(trim <= size());
        if (trim > 0) {
            rbegin += trim;
            if (qorientation) {
                qend -= trim;
            } else {
                qbegin += trim;
            }
            assert_valid(allow_empty);
            left_trim += trim;
        }
    }

    void trim_suffix(size_t trim, bool allow_empty = true) {
        assert_valid();
        assert(trim <= size());
        if (trim > 0) {
            rend -= trim;
            if (qorientation) {
                qbegin += trim;
            } else {
                qend -= trim;
            }
            assert_valid(allow_empty);
            right_trim += trim;
        }
    }

    std::string_view get_seq_from_target(std::string_view target) const {
        assert_valid();
        assert(rbegin + size() <= target.size());
        return std::string_view(target.data() + rbegin, size());
    }

    std::string_view get_seq_from_query(std::string_view query,
                                        std::string_view query_rc) const {
        assert_valid();
        assert(qorientation || qbegin + size() <= query.size());
        assert(!qorientation || static_cast<ssize_t>(query_rc.size()) >= qend);
        assert(!qorientation || query_rc.size() - qend + size() <= query_rc.size());
        return std::string_view(!qorientation
                                        ? query.data() + qbegin
                                        : query_rc.data() + query_rc.size() - qend,
                                 size());
    }

    bool check_equal(std::string_view target,
                     std::string_view query,
                     std::string_view query_rc) const {
        assert_valid();
        return get_seq_from_target(target) == get_seq_from_query(query, query_rc);
    }
};

inline bool operator<(const Ranges &a, const Ranges &b) {
    bool in_inv = a.qorientation && b.qorientation;

    return std::tie(a.rend, a.rbegin, !in_inv ? a.qend : b.qbegin, !in_inv ? a.qbegin : b.qend)
         < std::tie(b.rend, b.rbegin, !in_inv ? b.qend : a.qbegin, !in_inv ? b.qbegin : a.qend);
}

inline bool operator==(const Ranges &a, const Ranges &b) {
    return std::tie(a.rbegin, a.rend, a.rorientation, a.qbegin, a.qend, a.qorientation)
            == std::tie(b.rbegin, b.rend, b.rorientation, b.qbegin, b.qend, b.qorientation);
}

inline std::ostream& operator<<(std::ostream &out, const Ranges &a) {
    out << orientation_to_char(a.rorientation) << " "
        << a.rbegin + 1 << "-" << a.rend << "\t"
        << orientation_to_char(a.qorientation) << " "
        << a.qbegin + 1 << "-" << a.qend << "\t"
        << a.size() << "(" << a.left_trim << "," << a.right_trim << ")\t"
        << a.rend - a.qend;

    return out;
}

std::vector<Ranges> rc_ranges(const std::vector<Ranges>& mummer_ranges,
                              SOffset query_size);


struct ScoreModel {
    enum ModelType {
        EDIT_DISTANCE,
        GAP_LINEAR,
        GAP_AFFINE,
        GAP_2P_AFFINE
    };
    static const std::array<std::string, 4> model_type_str;

    static constexpr Penalty inf_p = std::numeric_limits<Penalty>::max() / 4;
    static constexpr Score ninf_s = std::numeric_limits<Score>::min() / 4;

    Score match_s;
    Score mismatch_s;
    Score gap_open_s;
    Score gap_ext_s;
    Score gap_open2_s;
    Score gap_ext2_s;
    Score inv_open_s;
    Score inv_ext_s;

    // the inversion extension score is divided by this. In practice, we implement
    // this by multiplying everything else
    Score inv_div_factor;

    ModelType model_type;

    Penalty mismatch_p;
    Penalty gap_open_p;
    Penalty gap_ext_p;
    Penalty gap_open2_p;
    Penalty gap_ext2_p;
    Penalty inv_open_p;
    Penalty inv_ext_p;

    SOffset gap_switch;
    Penalty penalty_diff;
    Penalty max_pen;

    ScoreModel(const ScoreModel&) = default;
    ScoreModel(ScoreModel&&) = default;
    ScoreModel& operator=(const ScoreModel&) = default;
    ScoreModel& operator=(ScoreModel&&) = default;

    ScoreModel()
          : match_s(2), mismatch_s(-2),
            gap_open_s(0), gap_ext_s(-3),
            gap_open2_s(0), gap_ext2_s(-3),
            inv_open_s(0), inv_ext_s(0), inv_div_factor(1),
            model_type(ModelType::EDIT_DISTANCE),
            mismatch_p(8),
            gap_open_p(0), gap_ext_p(8),
            gap_open2_p(0), gap_ext2_p(8),
            inv_open_p(0), inv_ext_p(0),
            gap_switch(max_offset),
            penalty_diff(inf_p),
            max_pen(8) {
        // edit distance by default
        // mismatch_p = (match_s - mismatch_s) * 2;
        // gap_ext_p = match_s - gap_ext_s * 2;

        // let mismatch_p = gap_ext_p = 8
        // then match_s - mismatch_s = 4, where match_s > 0 and mismatch_s < 0
        // so match_s = 2, mismatch_s = -2
        // so 8 = 2 - gap_ext_s * 2
        // so gap_ext_s = (2 - 8) / 2 = -3

        assert(match_s > 0);
        assert(mismatch_s < 0);
        assert(mismatch_p == static_cast<Penalty>(match_s - mismatch_s) * 2);
        assert(match_s > gap_ext_s * 2);
        assert(gap_ext_p == static_cast<Penalty>(match_s - gap_ext_s * 2));
    }

    ScoreModel(Score m, Score x, Score o1, Score e1, Score o2, Score e2, Score io, Score ie, Score idf = 1)
          : match_s(m), mismatch_s(x),
            gap_open_s(o1), gap_ext_s(e1),
            gap_open2_s(o2), gap_ext2_s(e2),
            inv_open_s(io), inv_ext_s(ie),
            inv_div_factor(idf),
            gap_switch(max_offset),
            penalty_diff(inf_p) {
        assert(match_s > 0);
        assert(mismatch_s < 0);
        assert(gap_open_s <= 0);
        assert(gap_ext_s < 0);
        assert(gap_open2_s <= 0);
        assert(gap_ext2_s < 0);

        assert(match_s > gap_ext_s * 2);

        assert(match_s > gap_open_s + gap_ext_s);
        assert(match_s > gap_open2_s + gap_ext2_s);

        assert(inv_open_s <= 0);
        assert(inv_ext_s <= 0);
        assert(match_s >= inv_open_s + inv_ext_s);

        if (gap_open_s != gap_open2_s || gap_ext_s != gap_ext2_s) {
            model_type = ModelType::GAP_2P_AFFINE;
        } else if (gap_open_s < 0) {
            model_type = ModelType::GAP_AFFINE;
        } else {
            model_type = ModelType::GAP_LINEAR;
        }

        if (gap_open2_s + gap_ext2_s > gap_open_s + gap_ext_s) {
            std::cerr << "WARNING: Swapping O1,E1 and O2,E2\n";
            std::swap(gap_open_s, gap_open2_s);
            std::swap(gap_ext_s, gap_ext2_s);
        }

        match_s *= inv_div_factor;
        mismatch_s *= inv_div_factor;
        gap_open_s *= inv_div_factor;
        gap_ext_s *= inv_div_factor;
        gap_open2_s *= inv_div_factor;
        gap_ext2_s *= inv_div_factor;
        inv_open_s *= inv_div_factor;

        // convert to equivalent penalties where match_p == 0
        // https://www.biorxiv.org/content/10.1101/2022.01.12.476087v1.full
        mismatch_p = (match_s - mismatch_s) * 2;
        assert(static_cast<Score>(mismatch_p) >= 0);

        gap_open_p = gap_open_s * -2;
        gap_ext_p = match_s - gap_ext_s * 2;
        assert(static_cast<Score>(gap_open_p) >= 0);
        assert(static_cast<Score>(gap_ext_p) > 0);

        gap_open2_p = gap_open2_s * -2;
        gap_ext2_p = match_s - gap_ext2_s * 2;
        assert(static_cast<Score>(gap_open2_p) >= 0);
        assert(static_cast<Score>(gap_ext2_p) > 0);

        // with inversion scoring
        // 2score + 2n_sN + 2i_sI = match_s*total_len - penalty + 2n_sN + 2i_sI
        // let 2n_s = -n_p and 2i_s = -i_p
        // then
        // 2score + 2n_sN + 2i_sI = match_s*total_len - penalty - n_pN - i_pI

        inv_open_p = inv_open_s * -2;
        inv_ext_p = inv_ext_s * -2;
        assert(static_cast<Score>(inv_open_p) >= 0);
        assert(static_cast<Score>(inv_ext_p) >= 0);

        if (gap_open_p < gap_open2_p && gap_ext_p > gap_ext2_p) {
            // there is a switchover
            // o1 + Le1 >= o2 + Le2
            // L(e1-e2) >= o2 - o1
            // L >= (o2 - o1) / (e1 - e2)
            gap_switch = std::ceil(static_cast<double>(gap_open2_p - gap_open_p) / (gap_ext_p - gap_ext2_p));
        }

        max_pen = std::max({ mismatch_p, gap_open_p + gap_ext_p, inv_open_p + inv_ext_p });
        if (gap_switch <= 1) {
            std::cerr << "WARNING: Swapping O1,E1 and O2,E2\n";
            gap_switch = max_offset;
            std::swap(gap_open_s, gap_open2_s);
            std::swap(gap_ext_s, gap_ext2_s);
            std::swap(gap_open_p, gap_open2_p);
            std::swap(gap_ext_p, gap_ext2_p);
        } else if (gap_switch < max_offset) {
            assert(gap_open2_p != gap_open_p || gap_ext2_p != gap_ext_p);
            penalty_diff = static_cast<Score>(gap_open2_p + gap_ext2_p * gap_switch)
                            - static_cast<Score>(gap_open_p + gap_ext_p * (gap_switch - 1));
            max_pen = std::max(max_pen, penalty_diff);
        }
    }

    ScoreModel(Score m, Score x, Score o, Score e, Score io, Score ie, Score idf = 1)
        : ScoreModel(m, x, o, e, o, e, io, ie, idf) {}

    const std::string& get_model_str() const {
        assert(model_type < model_type_str.size());
        return model_type_str[model_type];
    }

    Score penalty_to_score(Score p, SOffset qlen, SOffset rlen) const {
        Score init_penalty = match_s * (qlen + rlen) - p;
        assert(!(init_penalty % 2));
        return init_penalty / 2;
    }

    Score get_gap_score(SOffset gap_length) const {
        return gap_length > 0
            ? std::max(gap_open_s + gap_ext_s * gap_length,
                       gap_open2_s + gap_ext2_s * gap_length)
            : 0;
    }

    Penalty get_gap_penalty(SOffset gap_length) const {
        return gap_length > 0
            ? std::min(gap_open_p + gap_ext_p * gap_length,
                       gap_open2_p + gap_ext2_p * gap_length)
            : 0;
    }

    Penalty score_to_penalty(Score s, SOffset qlen, SOffset rlen) const {
        // s = (match_s * (qlen + rlen) - p) / 2
        // 2s = match_s * (qlen + rlen) - p
        // p = match_s * (qlen + rlen) - 2s
        Score init_penalty = match_s * (qlen + rlen) - s * 2;
        assert(init_penalty >= 0);
        return init_penalty;
    }
};

inline std::ostream& operator<<(std::ostream& out, const ScoreModel &model) {
    out << "Scores"
        << "\tModel: " << model.get_model_str()
        << "\tM: " << model.match_s / model.inv_div_factor
        << "\tX: " << model.mismatch_s / model.inv_div_factor
        << "\tO1: " << model.gap_open_s / model.inv_div_factor
        << ",E1: " << model.gap_ext_s / model.inv_div_factor
        << "\tO2: " << model.gap_open2_s / model.inv_div_factor
        << ",E2: " << model.gap_ext2_s / model.inv_div_factor
        << "\tIO: " << model.inv_open_s / model.inv_div_factor
        << ",IE: " << static_cast<double>(model.inv_ext_s) / model.inv_div_factor
        << "\nPenalties"
        << "\tM: " << 0 << "\tMM: " << model.mismatch_p
        << "\tO1: " << model.gap_open_p << "\tE1: " << model.gap_ext_p
        << "\tO2: " << model.gap_open2_p << "\tE2: " << model.gap_ext2_p
        << "\tIO: " << model.inv_open_p << "\tIE: " << model.inv_ext_p
        << "\tMaxPenaltyJump: " << model.max_pen
        << "\tGapModelSwitchLen: " << model.gap_switch
        << "\tGapModelSwitchPen: " << model.penalty_diff;
    return out;
}