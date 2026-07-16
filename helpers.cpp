#include "helpers.hpp"

#include <array>
#include <numeric>

#include <cctype>

std::array<char, 256> generate_rc_map() {
    std::array<char, 256> rc_map;
    std::iota(rc_map.begin(), rc_map.end(), 0);

    rc_map['A'] = 'T';
    rc_map['T'] = 'A';

    rc_map['G'] = 'C';
    rc_map['C'] = 'G';

    // A or G
    rc_map['R'] = 'Y';
    rc_map['Y'] = 'R';

    // G or T
    rc_map['K'] = 'M';
    rc_map['M'] = 'K';

    // C or G or T
    rc_map['B'] = 'V';
    rc_map['V'] = 'B';

    // A or G or T
    rc_map['D'] = 'H';
    rc_map['H'] = 'D';

    assert(rc_map['N'] == 'N');
    // G or C: S maps to S
    // A or T: W maps to W
    assert(rc_map['S'] == 'S');
    assert(rc_map['W'] == 'W');

    return rc_map;
}

static std::array<char, 256> rc_map = generate_rc_map();

std::string reverse_complement(std::string_view fw) {
    assert(std::all_of(fw.begin(), fw.end(), [](char c) { return c == std::toupper(c); }));
    std::string rc;
    rc.reserve(fw.size());
    std::transform(fw.rbegin(), fw.rend(), std::back_inserter(rc),
                   [&](char c) { return rc_map[c]; });
    return rc;
}