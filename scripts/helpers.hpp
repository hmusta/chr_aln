#pragma once

// Shared by syri_to_cigar and compute_edit_distance: FASTA reading, CIGAR
// operations, and the small parsing utilities both need.

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// ------------------------------------------------------------- parsing -----

inline std::vector<std::string_view> split(std::string_view line, char sep) {
    std::vector<std::string_view> fields;
    size_t begin = 0;
    while (true) {
        const size_t pos = line.find(sep, begin);
        if (pos == std::string_view::npos) {
            fields.push_back(line.substr(begin));
            break;
        }
        fields.push_back(line.substr(begin, pos - begin));
        begin = pos + 1;
    }
    return fields;
}

inline std::vector<std::string_view> split_ws(std::string_view line) {
    std::vector<std::string_view> fields;
    size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) ++i;
        if (i >= line.size()) break;
        const size_t begin = i;
        while (i < line.size() && line[i] != ' ' && line[i] != '\t') ++i;
        fields.push_back(line.substr(begin, i - begin));
    }
    return fields;
}

inline int64_t to_int(std::string_view s) {
    int64_t value = 0;
    bool negative = false;
    size_t i = 0;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) {
        negative = s[i] == '-';
        ++i;
    }
    if (i >= s.size())
        throw std::runtime_error("expected a number, got \"" + std::string(s) + "\"");
    for (; i < s.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(s[i])))
            throw std::runtime_error("expected a number, got \"" + std::string(s) + "\"");
        value = value * 10 + (s[i] - '0');
    }
    return negative ? -value : value;
}

inline void strip_eol(std::string *line) {
    while (!line->empty() && (line->back() == '\r' || line->back() == '\n'))
        line->pop_back();
}

// ------------------------------------------------------------ sequence -----

inline char complement(char c) {
    switch (c) {
        case 'A': return 'T';  case 'C': return 'G';
        case 'G': return 'C';  case 'T': return 'A';
        case 'R': return 'Y';  case 'Y': return 'R';
        case 'S': return 'S';  case 'W': return 'W';
        case 'K': return 'M';  case 'M': return 'K';
        case 'B': return 'V';  case 'V': return 'B';
        case 'D': return 'H';  case 'H': return 'D';
        default:  return c;    // N and anything unexpected map to themselves
    }
}

inline int64_t count_N(std::string_view window) {
    int64_t n = 0;
    for (const char c : window)
        n += (c == 'N');
    return n;
}

// Read a single-sequence FASTA, upper-casing as we go.
inline std::string read_fasta(const std::string &fname, std::string *header) {
    std::ifstream in(fname, std::ios::binary);
    if (!in)
        throw std::runtime_error("cannot open " + fname);

    std::ostringstream buffer;
    buffer << in.rdbuf();
    const std::string contents = buffer.str();

    const size_t pos = contents.find('\n');
    if (pos == std::string::npos)
        throw std::runtime_error(fname + ": no sequence after the header");

    std::string head = contents.substr(0, pos);
    if (!head.empty() && head[0] == '>')
        head.erase(0, 1);
    head = head.substr(0, head.find_first_of(" \t\r"));
    *header = head;

    std::string seq;
    seq.reserve(contents.size() - pos);
    for (size_t i = pos + 1; i < contents.size(); ++i) {
        const char c = contents[i];
        if (c != '\n' && c != '\r')
            seq.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return seq;
}

// --------------------------------------------------------------- CIGAR -----
//
// The convention both tools share:
//     =   identical bases, neither an N, consumes reference and query
//     X   differing bases, neither an N, consumes reference and query
//     M   a position where either side is an N, consumes reference and query
//     D   consumes reference only
//     I   consumes query only

struct Op {
    char op;
    int64_t length;
};

inline bool consumes_ref(char op) {
    return op == '=' || op == 'X' || op == 'M' || op == 'D';
}

inline bool consumes_qry(char op) {
    return op == '=' || op == 'X' || op == 'M' || op == 'I';
}

inline bool aligns_both(char op) {
    return op == '=' || op == 'X' || op == 'M';
}

// Append a run, merging it with the previous one when they share an operation.
inline void append_op(std::vector<Op> &ops, char op, int64_t length) {
    if (length <= 0)
        return;
    if (!ops.empty() && ops.back().op == op)
        ops.back().length += length;
    else
        ops.push_back({op, length});
}

// Split a CIGAR string into runs, keeping them exactly as written. Runs are
// deliberately not merged here: a length cutoff applied to an indel depends on
// how the aligner split it.
inline std::vector<Op> parse_cigar(std::string_view text) {
    std::vector<Op> ops;
    int64_t length = 0;
    bool seen_digit = false;
    for (const char c : text) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            length = length * 10 + (c - '0');
            seen_digit = true;
        } else {
            if (!seen_digit)
                throw std::runtime_error(std::string("CIGAR operation '") + c +
                                         "' has no length");
            ops.push_back({c, length});
            length = 0;
            seen_digit = false;
        }
    }
    if (seen_digit)
        throw std::runtime_error("CIGAR ends with a length and no operation");
    return ops;
}

inline std::pair<int64_t, int64_t> ops_span(const std::vector<Op> &ops) {
    int64_t r = 0, q = 0;
    for (const Op &o : ops) {
        if (consumes_ref(o.op)) r += o.length;
        if (consumes_qry(o.op)) q += o.length;
    }
    return {r, q};
}

// The CIGAR is the last whitespace-separated field of the first line.
inline std::string read_cigar(const std::string &fname) {
    std::ifstream in(fname);
    if (!in)
        throw std::runtime_error("cannot open " + fname);

    std::string line;
    if (!std::getline(in, line))
        return {};
    strip_eol(&line);

    const size_t pos = line.find_last_of(" \t");
    return pos == std::string::npos ? line : line.substr(pos + 1);
}
