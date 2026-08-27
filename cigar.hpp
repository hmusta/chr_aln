#pragma once

#include "helpers.hpp"

#include <functional>
#include <limits>
#include <ostream>
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

struct Cigar {
    using Storage = std::vector<std::pair<char, Offset>>;
    using iterator = typename Storage::iterator;
    using const_iterator = typename Storage::const_iterator;

    Cigar() = default;
    Cigar(const Cigar&) = default;
    Cigar(Cigar&&) = default;
    Cigar& operator=(const Cigar&) = default;
    Cigar& operator=(Cigar&&) = default;

    Cigar(char c, Offset num) : data_(1, std::make_pair(c, num)) {}
    Cigar(const std::string &cigar_str) {
        Offset num = 0;
        for (unsigned char c : cigar_str) {
            if (std::isdigit(c)) {
                num = num * 10 + c - '0';
            } else {
                assert(num);
                data_.emplace_back(c, num);
                num = 0;
            }
        }
    }

    bool empty() const { return data_.empty(); }
    void clear() { data_.clear(); }

    iterator begin() { return data_.begin(); }
    const_iterator begin() const { return data_.begin(); }
    const_iterator cbegin() const { return data_.cbegin(); }

    iterator end() { return data_.end(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cend() const { return data_.cend(); }

    void push(char c, Offset num, bool merge = false) {
        switch (c) {
            case EQ_OP:
            case MATCH_OP:
            case NEQ_OP:
            case TARGET_CONSUME_OP:
            case QUERY_CONSUME_OP:
            case INV_OP: {
                if (merge && data_.size() && data_.back().first == c) {
                    data_.back().second += num;
                } else {
                    data_.emplace_back(c, num);
                }
            } break;
            default: {
                std::cerr << "\n\nInvalid character pushed to in CIGAR\n"
                        << c << "\n\n"
                        << std::endl;
                assert(false);
            }
        }
    }

    template <typename It>
    void insert(iterator it, It begin, It end, bool merge = false) {
        if (merge) {
            if (it != data_.end() && begin != end && (end - 1)->first == it->first) {
                it->second += (end - 1)->second;
                --end;
            }

            if (it != data_.begin() && begin != end && begin->first == (it - 1)->first) {
                (it - 1)->second += begin->second;
                ++begin;
            }
        }

        data_.insert(it, begin, end);
    }

    bool operator==(const Cigar &b) const { return data_ == b.data_; }

    Storage data_;
};

inline std::ostream& operator<<(std::ostream& out, const Cigar &cigar) {
    for (const auto &[c, num] : cigar) {
        out << num << c;
    }
    return out;
}

inline void cigar_caller(const Cigar& cigar,
                         const std::function<void(char, Offset, Offset, Offset)>& callback,
                         Offset r_len = std::numeric_limits<Offset>::max(),
                         Offset q_len = std::numeric_limits<Offset>::max(),
                         const std::function<bool()>& terminate = []() { return false; }) {
    Offset r_pos = 0;
    Offset q_pos = 0;

    for (auto [c, num] : cigar) {
        if (terminate())
            return;

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
    }

    assert(r_len == std::numeric_limits<Offset>::max() || r_len == r_pos);
    assert(q_len == std::numeric_limits<Offset>::max() || q_len == q_pos);

    #ifdef NDEBUG
    std::ignore = r_len;
    std::ignore = q_len;
    #endif
}

inline std::pair<int64_t, int64_t> count_identities_and_matches(const Cigar& cigar,
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

inline Score score_cigar(const Cigar& cigar,
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

inline size_t cigar_edits(const Cigar& cigar,
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

inline std::pair<Offset, Offset>
cigar_get_target_pos(const Cigar& cigar,
                     size_t final_query_pos,
                     Offset r_len,
                     Offset q_len) {
    assert(final_query_pos <= q_len);
    if (final_query_pos == q_len)
        return std::make_pair(r_len, r_len);

    bool terminate = false;
    Offset target_pos_start = 0;
    Offset target_pos_end = 0;
    cigar_caller(cigar, [&](char c, Offset num, Offset r_pos, Offset q_pos) {
        assert(q_pos <= final_query_pos);
        assert(target_pos_end == r_pos);
        Offset old_num = num;

        switch (c) {
            case EQ_OP:
            case MATCH_OP:
            case NEQ_OP: {
                if (q_pos >= final_query_pos) {
                    terminate = true;
                    return;
                }

                size_t next_query_pos = std::min<Offset>(final_query_pos, q_pos + num);
                num = next_query_pos - q_pos;
                target_pos_start = r_pos + num;
                target_pos_end = target_pos_start;
                terminate = (next_query_pos == final_query_pos && old_num > num);
            } break;
            case QUERY_CONSUME_OP: {
                if (q_pos >= final_query_pos) {
                    terminate = true;
                    return;
                }

                terminate = (q_pos + num > final_query_pos);
            } break;
            case TARGET_CONSUME_OP: {
                if (q_pos >= final_query_pos) {
                    target_pos_end += num;
                } else {
                    target_pos_start += num;
                    target_pos_end = target_pos_start;
                }
            } break;
        }
    }, r_len, q_len, [&]() { return terminate; });

    return std::make_pair(target_pos_start, target_pos_end);
}

inline std::pair<Offset, Offset>
cigar_get_query_pos(const Cigar& cigar,
                     size_t final_target_pos,
                     Offset r_len,
                     Offset q_len) {
    assert(final_target_pos <= r_len);
    if (final_target_pos == r_len)
        return std::make_pair(q_len, q_len);

    bool terminate = false;
    Offset query_pos_start = 0;
    Offset query_pos_end = 0;
    cigar_caller(cigar, [&](char c, Offset num, Offset r_pos, Offset q_pos) {
        assert(r_pos <= final_target_pos);
        assert(query_pos_end == q_pos);
        Offset old_num = num;

        switch (c) {
            case EQ_OP:
            case MATCH_OP:
            case NEQ_OP: {
                if (r_pos >= final_target_pos) {
                    terminate = true;
                    return;
                }

                size_t next_target_pos = std::min<Offset>(final_target_pos, r_pos + num);
                num = next_target_pos - r_pos;
                query_pos_start = q_pos + num;
                query_pos_end = query_pos_start;
                terminate = (next_target_pos == final_target_pos && old_num > num);
            } break;
            case TARGET_CONSUME_OP: {
                if (r_pos >= final_target_pos) {
                    terminate = true;
                    return;
                }

                terminate = (r_pos + num > final_target_pos);
            } break;
            case QUERY_CONSUME_OP: {
                if (r_pos >= final_target_pos) {
                    query_pos_end += num;
                } else {
                    query_pos_start += num;
                    query_pos_end = query_pos_start;
                }
            } break;
        }
    }, r_len, q_len, [&]() { return terminate; });

    return std::make_pair(query_pos_start, query_pos_end);
}

inline std::pair<Cigar, Cigar>
cigar_split(const Cigar& cigar,
            size_t final_target_pos,
            size_t final_query_pos,
            Offset r_len,
            Offset q_len) {
    assert(final_target_pos <= r_len);
    assert(final_query_pos <= q_len);
    if (final_target_pos == r_len && final_query_pos == q_len)
        return std::make_pair(cigar, Cigar());

    size_t target_pos = 0;
    size_t query_pos = 0;
    auto terminate = [&]() {
        assert(target_pos <= r_len);
        assert(query_pos <= q_len);
        return target_pos == final_target_pos && query_pos == final_query_pos;
    };

    Cigar prefix;
    Cigar suffix;
    auto it = cigar.begin();

    cigar_caller(cigar, [&](char c, Offset num, Offset r_pos, Offset q_pos) {
        assert(it != cigar.end());
        assert(q_pos <= final_query_pos);
        assert(r_pos <= final_target_pos);
        assert(target_pos == r_pos);
        assert(query_pos == q_pos);

        Offset old_num = num;

        switch (c) {
            case EQ_OP:
            case MATCH_OP:
            case NEQ_OP: {
                target_pos = std::min<size_t>(final_target_pos, r_pos + num);
                query_pos = std::min<size_t>(final_query_pos, q_pos + num);
                num = std::min(target_pos - r_pos, query_pos - q_pos);
                assert(num > 0);
            } break;
            case QUERY_CONSUME_OP: {
                query_pos = std::min<size_t>(final_query_pos, q_pos + num);
                num = query_pos - q_pos;
                assert(num > 0);
            } break;
            case TARGET_CONSUME_OP: {
                target_pos = std::min<size_t>(final_target_pos, r_pos + num);
                num = target_pos - r_pos;
                assert(num > 0);
            } break;
        }

        prefix.push(c, num);

        if (old_num > num)
            suffix = Cigar(c, old_num - num);

        ++it;
    }, r_len, q_len, terminate);

    suffix.insert(suffix.end(), it, cigar.end());

    return std::make_pair(std::move(prefix), std::move(suffix));
}

inline Cigar cigar_fix_n(const Cigar& cigar_in,
                           std::string_view target,
                           std::string_view query) {
    Cigar cigar;

    auto push_op = [&](char c, int64_t num) {
        assert(c == MATCH_OP || c == NEQ_OP || c == EQ_OP);
        cigar.push(c, num);
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
                cigar.push(c, num);
                last_op = c;
            } break;
        }
    }, target.size(), query.size());

    if (snum > 0)
        push_op(last_op, snum);

    return cigar;
}
