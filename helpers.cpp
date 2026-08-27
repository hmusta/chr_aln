#include "helpers.hpp"

#include <array>
#include <numeric>

#include <cctype>

const std::array<std::string, 4> ScoreModel::model_type_str {
        "Edit distance", "Gap linear", "Gap affine", "2-piece gap affine"
    };

static_assert((1llu << (sizeof(unsigned char) * 8)) == 256);
static std::array<unsigned char, 256> RC_MAP = []() {
    std::array<unsigned char, 256> rc_map;

    // by default, the identity map
    std::iota(rc_map.begin(), rc_map.end(), 0);

    rc_map['A'] = 'T';
    rc_map['a'] = 't';
    rc_map['T'] = 'A';
    rc_map['t'] = 'a';

    rc_map['G'] = 'C';
    rc_map['g'] = 'c';
    rc_map['C'] = 'G';
    rc_map['c'] = 'g';

    // A or G <-> T or C
    rc_map['R'] = 'Y';
    rc_map['r'] = 'y';
    rc_map['Y'] = 'R';
    rc_map['y'] = 'r';

    // G or T <-> C or A
    rc_map['K'] = 'M';
    rc_map['k'] = 'm';
    rc_map['M'] = 'K';
    rc_map['m'] = 'k';

    // C or G or T <-> G or C or A
    rc_map['B'] = 'V';
    rc_map['b'] = 'v';
    rc_map['V'] = 'B';
    rc_map['v'] = 'b';

    // A or G or T <-> T or C or A
    rc_map['D'] = 'H';
    rc_map['d'] = 'h';
    rc_map['H'] = 'D';
    rc_map['h'] = 'd';

    // all others should be the identity, but double-check the remaining characters
    assert(rc_map['N'] == 'N');
    assert(rc_map['n'] == 'n');

    // G or C <-> C or G
    assert(rc_map['S'] == 'S');
    assert(rc_map['s'] == 's');

    // A or T <-> T or A
    assert(rc_map['W'] == 'W');
    assert(rc_map['w'] == 'w');

    return rc_map;
}();

std::string reverse_complement(std::string_view fw) {
    assert(std::all_of(fw.begin(), fw.end(), [](unsigned char c) { return c == std::toupper(c); }));

    std::string rc;
    rc.reserve(fw.size());
    std::transform(fw.rbegin(), fw.rend(), std::back_inserter(rc),
                   [](unsigned char c) { return RC_MAP[c]; });
    return rc;
}

bool is_reverse_complement(std::string_view fw, std::string_view rc) {
    return std::equal(fw.rbegin(), fw.rend(), rc.begin(), rc.end(), [](unsigned char f, unsigned char r) {
        assert(f == std::toupper(f));
        assert(r == std::toupper(r));
        return f == RC_MAP[r];
    });
}

std::vector<Ranges> rc_ranges(const std::vector<Ranges>& mummer_ranges,
                              SOffset query_size) {
    std::vector<Ranges> rc_mummer_ranges;
    rc_mummer_ranges.reserve(mummer_ranges.size());
    for (const auto& [rbegin, rend, rrc, qbegin, qend, qrc, left_trim, right_trim] : mummer_ranges) {
        assert(query_size >= qbegin);
        assert(query_size >= qend);
        rc_mummer_ranges.emplace_back(rbegin, rend, rrc, query_size - qend,
                                      query_size - qbegin, !qrc, left_trim, right_trim);
    }

    return rc_mummer_ranges;
}