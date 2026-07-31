#pragma once

#include "helpers.hpp"

#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <utility>

#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstring>

// canonical CIGAR operations
// static const char *CIGAR_OPS = "MIDNSHP=X";

static constexpr char MATCH_OP = 'M'; // consumes both
static constexpr char QUERY_CONSUME_OP = 'I';
static constexpr char TARGET_CONSUME_OP = 'D';
static constexpr char TARGET_CLIP_OP = 'N';
static constexpr char QUERY_CLIP_OP = 'S';
static constexpr char QUERY_HARDCLIP_OP = 'H';
static constexpr char PAD_OP = 'P'; // consumes none
static constexpr char EQ_OP = '='; // consumes both
static constexpr char NEQ_OP = 'X'; // consumes both
static constexpr char INV_OP = 'i'; // consumes none, but query should have space for it
static constexpr char N_OP = '\0';

inline bool consumes_target(char c) {
    switch (c) {
        case MATCH_OP:
        case EQ_OP:
        case NEQ_OP:
        case TARGET_CONSUME_OP: { return true; } break;
        case QUERY_CONSUME_OP:
        case INV_OP: { return false; } break;
    }

    std::cerr << "\n\nInvalid CIGAR char\n"
            << c << "\n\n"
            << std::endl;
    assert(false);
    return false;
}

inline bool consumes_query(char c) {
    switch (c) {
        case MATCH_OP:
        case EQ_OP:
        case NEQ_OP:
        case QUERY_CONSUME_OP: { return true; } break;
        case TARGET_CONSUME_OP:
        case INV_OP: { return false; } break;
    }

    std::cerr << "\n\nInvalid CIGAR char\n"
            << c << "\n\n"
            << std::endl;
    assert(false);
    return false;
}

inline void cigar_caller(const std::string& cigar,
                         const std::function<void(char, Offset, Offset, Offset)>& callback,
                         Offset r_len = std::numeric_limits<Offset>::max(),
                         Offset q_len = std::numeric_limits<Offset>::max(),
                         const std::function<bool()>& terminate = []() { return false; }) {
    assert(cigar.empty() || cigar.size() >= 2);
    assert(cigar.empty() || std::isdigit(cigar[0]));
    assert(cigar.empty() || !std::isdigit(cigar.back()));

    Offset num = 0;
    Offset r_pos = 0;
    Offset q_pos = 0;

    for (char c : cigar) {
        if (terminate())
            return;

        if (std::isdigit(c)) {
            num = num * 10 + c - '0';
        } else {
            assert(num);
            if (!terminate()) {
                switch (c) {
                    case EQ_OP:
                    case MATCH_OP:
                    case NEQ_OP: {
                        assert(r_pos + num <= r_len);
                        assert(q_pos + num <= q_len);
                        callback(c, num, r_pos, q_pos);
                        r_pos += num;
                        q_pos += num;
                    } break;
                    case TARGET_CONSUME_OP: {
                        assert(r_pos + num <= r_len);
                        callback(c, num, r_pos, q_pos);
                        r_pos += num;
                    } break;
                    case QUERY_CONSUME_OP: {
                        assert(q_pos + num <= q_len);
                        callback(c, num, r_pos, q_pos);
                        q_pos += num;
                    } break;
                    case INV_OP: {
                        assert(q_pos + num <= q_len);
                        callback(c, num, r_pos, q_pos);
                    } break;
                    default: {
                        std::cerr << "\n\nInvalid character in CIGAR\n"
                                << cigar << "\n\n"
                                << std::endl;
                        assert(false);
                    }
                }
            } else {
                return;
            }
            num = 0;
        }
    }

    assert(num == 0);
    assert(r_len == std::numeric_limits<Offset>::max() || r_len == r_pos);
    assert(q_len == std::numeric_limits<Offset>::max() || q_len == q_pos);

    #ifdef NDEBUG
    std::ignore = r_len;
    std::ignore = q_len;
    #endif
}

inline std::pair<int64_t, int64_t> count_identities_and_matches(const std::string& cigar,
                                                                Offset r_len = std::numeric_limits<Offset>::max(),
                                                                Offset q_len = std::numeric_limits<Offset>::max()) {
    int64_t identities = 0;
    int64_t matches = 0;
    cigar_caller(cigar, [&](char c, Offset num, Offset /* r_pos */, Offset /* q_pos */) {
        switch (c) {
            case EQ_OP: identities += num; [[fallthrough]];
            case NEQ_OP:
            case MATCH_OP: { matches += num; } break;
        }
    }, r_len, q_len);

    return std::make_pair(identities, matches);
}

inline Score score_cigar(const std::string& cigar,
                         const SeqPair& view_pair,
                         const ScoreModel &score_model,
                         bool penalty = false) {
    auto [query, target] = view_pair;
    auto it_q = query.begin();
    auto it_r = target.begin();
    Score score = 0;
    cigar_caller(cigar, [&](char c, Offset num, Offset r_pos, Offset q_pos) {
        assert(it_q == query.begin() + q_pos);
        assert(it_r == target.begin() + r_pos);
        #ifdef NDEBUG
        std::ignore = q_pos;
        std::ignore = r_pos;
        #endif
        switch (c) {
            case MATCH_OP: {
                if (!penalty)
                    score += score_model.match_s * num;

                assert(std::equal(it_r, it_r + num, it_q, [](char a, char b) { return a == 'N' || b == 'N'; }));
                it_q += num;
                it_r += num;
            } break;
            case EQ_OP: {
                if (!penalty)
                    score += score_model.match_s * num;
                assert(std::equal(it_q, it_q + num, it_r));
                it_q += num;
                it_r += num;
            } break;
            case NEQ_OP: {
                score += (!penalty ? score_model.mismatch_s : score_model.mismatch_p) * num;
                assert(std::equal(it_q, it_q + num, it_r, [](char a, char b) { return a != b; }));
                it_q += num;
                it_r += num;
            } break;
            case TARGET_CONSUME_OP: {
                score += !penalty ? score_model.get_gap_score(num) : score_model.get_gap_penalty(num);
                it_r += num;
            } break;
            case QUERY_CONSUME_OP: {
                score += !penalty ? score_model.get_gap_score(num) : score_model.get_gap_penalty(num);
                it_q += num;
            } break;
        }
    }, target.size(), query.size());

    return score;
}

inline size_t cigar_edits(const std::string& cigar,
                          Offset r_len = std::numeric_limits<Offset>::max(),
                          Offset q_len = std::numeric_limits<Offset>::max()) {
    size_t edits = 0;
    cigar_caller(cigar, [&](char c, Offset num, Offset /* r_pos */, Offset /* q_pos */) {
        switch (c) {
            case MATCH_OP:
            case EQ_OP: break;
            case NEQ_OP:
            case TARGET_CONSUME_OP:
            case QUERY_CONSUME_OP: { edits += num; } break;
        }
    }, r_len, q_len);

    return edits;
}

inline size_t cigar_get_target_pos(const std::string& cigar,
                                   size_t final_query_pos,
                                   Offset r_len,
                                   Offset q_len) {
    assert(final_query_pos <= q_len);
    bool terminate = (final_query_pos == q_len);
    size_t target_pos = 0;
    cigar_caller(cigar, [&](char c, Offset num, Offset r_pos, Offset q_pos) {
        assert(q_pos <= final_query_pos);
        assert(target_pos == r_pos);

        switch (c) {
            case EQ_OP:
            case MATCH_OP:
            case NEQ_OP: {
                size_t next_query_pos = std::min<size_t>(final_query_pos, q_pos + num);
                num = next_query_pos - q_pos;
                target_pos = r_pos + num;
                terminate = (next_query_pos == final_query_pos);
            } break;
            case QUERY_CONSUME_OP: {
                terminate = (q_pos + num >= final_query_pos);
            } break;
            case TARGET_CONSUME_OP: { target_pos += num; } break;
        }
    }, r_len, q_len, [&]() { return terminate; });

    return target_pos;
}

inline size_t cigar_get_query_pos(const std::string& cigar,
                                  size_t final_target_pos,
                                  Offset r_len,
                                  Offset q_len) {
    assert(final_target_pos <= r_len);
    bool terminate = (final_target_pos == r_len);
    size_t query_pos = 0;
    cigar_caller(cigar, [&](char c, Offset num, Offset r_pos, Offset q_pos) {
        assert(r_pos <= final_target_pos);
        assert(query_pos == q_pos);

        switch (c) {
            case EQ_OP:
            case MATCH_OP:
            case NEQ_OP: {
                size_t next_target_pos = std::min<size_t>(final_target_pos, r_pos + num);
                num = next_target_pos - r_pos;
                query_pos = q_pos + num;
                terminate = (next_target_pos == final_target_pos);
            } break;
            case QUERY_CONSUME_OP: { query_pos += num; } break;
            case TARGET_CONSUME_OP: {
                terminate = (r_pos + num >= final_target_pos);
            } break;
        }
    }, r_len, q_len, [&]() { return terminate; });

    return query_pos;
}

inline std::string cigar_fix_n(const std::string& cigar_in,
                               std::string_view target,
                               std::string_view query) {
    std::string cigar;

    auto push_op = [&](char c, int64_t num) {
        assert(c == MATCH_OP || c == NEQ_OP || c == EQ_OP);
        cigar += std::to_string(num) + c;
    };

    Offset snum = 0;
    char last_op = N_OP;
    cigar_caller(cigar_in, [&](char c, Offset num, Offset r_pos, Offset q_pos) {
        switch (c) {
            case MATCH_OP:
            case EQ_OP:
            case NEQ_OP: {
                std::string_view target_w = target.substr(r_pos, num);
                std::string_view query_w = query.substr(q_pos, num);
                for (size_t i = 0; i < target_w.size(); ++i) {
                    char op;
                    if (target_w[i] != 'N' && query_w[i] != 'N') {
                        op = (target_w[i] == query_w[i]) ? EQ_OP : NEQ_OP;
                    } else {
                        op = MATCH_OP;
                    }
                    if (op == last_op) {
                        ++snum;
                    } else {
                        if (snum > 0)
                            push_op(last_op, snum);

                        snum = 1;
                    }
                    last_op = op;
                }
            } break;
            case TARGET_CONSUME_OP:
            case QUERY_CONSUME_OP:
            case INV_OP: {
                if (snum > 0) {
                    push_op(last_op, snum);
                    snum = 0;
                    last_op = N_OP;
                }
                cigar += std::to_string(num) + c;
                last_op = c;
            } break;
        }
    }, target.size(), query.size());

    if (snum > 0)
        push_op(last_op, snum);

    return cigar;
}
