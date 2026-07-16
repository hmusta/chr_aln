#pragma once

#include "helpers.hpp"

#include <functional>
#include <string>
#include <string_view>
#include <utility>

#include <cassert>
#include <cctype>
#include <cstdint>

inline void cigar_caller(const std::string& cigar,
                  const std::function<void(char, int64_t)>& callback,
                  const std::function<bool()>& terminate = []() { return false; }) {
    int64_t num = 0;
    for (char c : cigar) {
        if (terminate())
            return;

        if (std::isdigit(c)) {
            num = num * 10 + c - '0';
        } else {
            assert(num);
            if (!terminate()) {
                callback(c, num);
            } else {
                return;
            }
            num = 0;
        }
    }

    assert(num == 0);
}

inline std::pair<int64_t, int64_t> count_identities_and_matches(const std::string& cigar) {
    int64_t identities = 0;
    int64_t matches = 0;
    cigar_caller(cigar, [&](char c, int64_t num) {
        switch (c) {
            case '=': identities += num; [[fallthrough]];
            case 'X':
            case 'M': { matches += num; } break;
        }
    });

    return std::make_pair(identities, matches);
}

inline bool check_cigar_seq_lengths(const std::string& cigar,
                                    SOffset r_consumed,
                                    SOffset q_consumed) {
    SOffset r_consumed_test = 0;
    SOffset q_consumed_test = 0;
    cigar_caller(cigar, [&](char c, int64_t num) {
        switch (c) {
            case 'X':
            case 'M':
            case '=': {
                assert(r_consumed_test + num <= r_consumed);
                assert(q_consumed_test + num <= q_consumed);
                r_consumed_test += num;
                q_consumed_test += num;
            } break;
            case 'D': {
                assert(q_consumed_test + num <= q_consumed);
                q_consumed_test += num;
            } break;
            case 'I': {
                assert(r_consumed_test + num <= r_consumed);
                r_consumed_test += num;
            } break;
            default: {
                std::cerr << "\n\nInvalid character in CIGAR\n"
                          << cigar << "\n\n"
                          << std::endl;
                assert(false);
            }
        }
    });

    assert(r_consumed_test == r_consumed);
    assert(q_consumed_test == q_consumed);
    return r_consumed_test == r_consumed && q_consumed_test == q_consumed;
}

inline Score score_cigar(const std::string& cigar,
                         const SeqPair& view_pair,
                         const ScoreModel &score_model,
                         bool penalty = false) {
    auto [query, target] = view_pair;
    auto it_q = query.begin();
    auto it_r = target.begin();
    Score score = 0;
    cigar_caller(cigar, [&](char c, int64_t num) {
        switch (c) {
            case 'M': {
                if (!penalty)
                    score += score_model.match_s * num;
                assert(it_q + num <= query.end());
                assert(it_r + num <= target.end());
                assert(std::equal(it_r, it_r + num, it_q, [](char a, char b) { return a == 'N' || b == 'N'; }));
                it_q += num;
                it_r += num;
            } break;
            case '=': {
                if (!penalty)
                    score += score_model.match_s * num;
                assert(it_q + num <= query.end());
                assert(it_r + num <= target.end());
                assert(std::equal(it_q, it_q + num, it_r));
                it_q += num;
                it_r += num;
            } break;
            case 'X': {
                score += (!penalty ? score_model.mismatch_s : score_model.mismatch_p) * num;
                assert(it_q + num <= query.end());
                assert(it_r + num <= target.end());
                #ifndef NDEBUG
                for (int64_t i = 0; i < num; ++i) {
                    assert(*it_q != *it_r);
                    ++it_q;
                    ++it_r;
                }
                #endif
            } break;
            case 'I': {
                score += !penalty ? score_model.get_gap_score(num) : score_model.get_gap_penalty(num);
                it_r += num;
                assert(it_r <= target.end());
            } break;
            case 'D': {
                score += !penalty ? score_model.get_gap_score(num) : score_model.get_gap_penalty(num);
                it_q += num;
                assert(it_q <= query.end());
            } break;
            default: {
                std::cerr << "\n\nInvalid character in CIGAR\n"
                            << cigar << "\n\n"
                            << std::endl;
                assert(false);
            }
        }
    });

    assert(it_r == target.end());
    assert(it_q == query.end());

    return score;
}

inline size_t cigar_edits(const std::string& cigar) {
    size_t edits = 0;
    cigar_caller(cigar, [&](char c, int64_t num) {
        switch (c) {
            case 'M':
            case '=': break;
            case 'X':
            case 'I':
            case 'D': { edits += num; } break;
            default: {
                std::cerr << "\n\nInvalid character in CIGAR\n"
                            << cigar << "\n\n"
                            << std::endl;
                assert(false);
            }
        }
    });

    return edits;
}

inline size_t cigar_get_target_pos(const std::string& cigar, size_t final_query_pos) {
    size_t target_pos = 0;
    size_t query_pos = 0;
    cigar_caller(cigar, [&](char c, int64_t num) {
        assert(query_pos <= final_query_pos);
        switch (c) {
            case '=':
            case 'M':
            case 'X': {
                size_t next_query_pos = std::min<size_t>(final_query_pos, query_pos + num);
                num = next_query_pos - query_pos;
                query_pos = next_query_pos;
                target_pos += num;
            } break;
            case 'D': {
                query_pos = std::min<size_t>(final_query_pos, query_pos + num);
            } break;
            case 'I': { target_pos += num; } break;
            default: {
                std::cerr << "\n\nInvalid character in CIGAR\n"
                            << cigar << "\n\n"
                            << std::endl;
                assert(false);
            }
        }
    }, [&]() { return query_pos >= final_query_pos; });

    return target_pos;
}

inline size_t cigar_get_query_pos(const std::string& cigar, size_t final_target_pos) {
    size_t target_pos = 0;
    size_t query_pos = 0;
    cigar_caller(cigar, [&](char c, int64_t num) {
        assert(target_pos <= final_target_pos);
        switch (c) {
            case '=':
            case 'M':
            case 'X': {
                size_t next_target_pos = std::min<size_t>(final_target_pos, target_pos + num);
                num = next_target_pos - target_pos;
                target_pos = next_target_pos;
                query_pos += num;
            } break;
            case 'D': { query_pos += num; } break;
            case 'I': {
                target_pos = std::min<size_t>(final_target_pos, target_pos + num);
            } break;
            default: {
                std::cerr << "\n\nInvalid character in CIGAR\n"
                            << cigar << "\n\n"
                            << std::endl;
                assert(false);
            }
        }
    }, [&]() { return target_pos >= final_target_pos; });

    return query_pos;
}

inline std::string cigar_fix_n(const std::string& cigar_in,
                               std::string_view target,
                               std::string_view query) {
    std::string cigar;
    Offset r_consumed = 0;
    Offset q_consumed = 0;

    auto push_op = [&](char c, int64_t num) {
        assert(c == 'M' || c == 'X' || c == '=');
        cigar += std::to_string(num) + c;
    };

    Offset snum = 0;
    char last_op = 'S';
    cigar_caller(cigar_in, [&](char op, int64_t num) {
        switch (op) {
            case 'M':
            case '=':
            case 'X': {
                assert(r_consumed + num <= target.size());
                assert(q_consumed + num <= query.size());
                std::string_view target_w = target.substr(r_consumed, num);
                std::string_view query_w = query.substr(q_consumed, num);
                for (size_t i = 0; i < target_w.size(); ++i) {
                    char op;
                    if (target_w[i] != 'N' && query_w[i] != 'N') {
                        op = (target_w[i] == query_w[i]) ? '=' : 'X';
                    } else {
                        op = 'M';
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
                r_consumed += num;
                q_consumed += num;
            } break;
            case 'I': {
                if (snum > 0) {
                    push_op(last_op, snum);
                    snum = 0;
                    last_op = 'S';
                }
                assert(r_consumed + num <= target.size());
                cigar += std::to_string(num) + 'I';
                r_consumed += num;
                last_op = 'I';
            } break;
            case 'D': {
                if (snum > 0) {
                    push_op(last_op, snum);
                    snum = 0;
                    last_op = 'S';
                }
                assert(q_consumed + num <= query.size());
                cigar += std::to_string(num) + 'D';
                q_consumed += num;
                last_op = 'D';
            } break;
        }
    });

    if (snum > 0)
        push_op(last_op, snum);

    assert(r_consumed == target.size());
    assert(q_consumed == query.size());

    assert(target.find('N') != std::string_view::npos
            || query.find('N') != std::string_view::npos
            || cigar == cigar_in);

    return cigar;
}