#pragma once

#include "helpers.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <cstdint>

inline std::pair<std::string, std::string> read_fasta(std::istream &fin) {
    std::string fasta;
    std::string line;
    std::string header_full;
    std::getline(fin, header_full);
    std::string header;
    {
        std::istringstream sin(header_full);
        sin >> header;
    }
    header = header.substr(1);
    while (fin >> line) {
        if (line[0] == '>')
            break;

        fasta += line;
    }

    std::transform(fasta.begin(), fasta.end(), fasta.begin(), sanitize_nuc);

    return std::make_pair(std::move(fasta), std::move(header));
}

inline std::vector<Ranges> read_mummer(std::istream &fin,
                                       std::string &target,
                                       std::string &theader,
                                       std::string &query,
                                       std::string &qheader,
                                       std::string &query_rc) {
    std::vector<Ranges> output;
    std::string line;

    bool qrc = false;
    while (std::getline(fin, line)) {
        if (line[0] == '>') {
            auto rev_idx = line.rfind("Reverse");
            bool is_rev = (rev_idx != std::string::npos);
            assert(rev_idx == std::string::npos || rev_idx > 3);
            std::string_view query_header(line.c_str() + 2,
                                          !is_rev ? line.size() - 2 : rev_idx - 3);
            assert(query_header == theader || query_header == qheader);
            if (is_rev)
                qrc = true;

            if (!is_rev && query_header == theader)
                throw std::runtime_error("ERROR: reference and query in the wrong order");

            continue;
        }

        auto add_mum = [&](SOffset rbegin, SOffset qbegin, SOffset len) {
            assert(rbegin > 0);
            assert(qbegin > 0);
            assert(target.size() >= static_cast<size_t>(rbegin));
            assert(query.size() >= static_cast<size_t>(qbegin));
            --rbegin;
            --qbegin;

            assert(static_cast<size_t>(rbegin) < target.size());
            assert(static_cast<size_t>(rbegin + len) <= target.size());
            assert(static_cast<size_t>(qbegin) < query.size());
            assert(static_cast<size_t>(qbegin + len) <= query.size());

            if (qrc) {
                // match is [query.size() - qbegin - 1, query.size() - qbegin - 1 + len) in query_rc
                // so [qbegin - len + 1, qbegin + 1) in query
                assert(qbegin + 1 >= len);
                qbegin -= len - 1;
            }

            output.emplace_back(rbegin, rbegin + len, false, qbegin, qbegin + len, qrc);
            assert(output.back().check_equal(target, query, query_rc));
        };

        std::istringstream sin(line);
        Offset rbegin;
        Offset qbegin;
        Offset len;

        sin >> rbegin >> qbegin >> len;

        add_mum(rbegin, qbegin, len);
    }

    return output;
}