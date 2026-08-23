// Summarise a global CIGAR alignment against the two sequences it aligns.
//
// Usage:
//     compute_edit_distance CHECK_N CHECK_EQ (REF.fa QRY.fa CIGAR INV.bed)...
//
// CHECK_N and CHECK_EQ are 0/1. Each group of four arguments describes one
// chromosome; INV.bed may be absent, in which case it is skipped.
//
// The CIGAR uses = and X only for positions where neither sequence has an N,
// and M for positions aligned against an N on either side. All three consume
// reference and query alike; D and I consume one side each.
//
// With both checks clear the operations are taken at their word: = is a match
// between two non-N bases, X a mismatch between two non-N bases, and M a
// position against an N. Nothing is read back from the sequences except the
// N counts inside M, D and I runs, which are needed to shorten the lengths.
//
// With either check set nothing is assumed and every aligned run is measured
// against the sequences instead: = must be identical throughout (CHECK_EQ), X
// must differ throughout (CHECK_EQ), and M must really cover an N on one side
// or the other. The run is then scored by what its bases actually are, so a
// mislabelled one is counted where it belongs rather than where its label
// would have put it. CHECK_N is what a CIGAR with no M convention needs, since
// there the Ns are hidden inside = and X runs and can only be found by
// looking.
//
// A C++ rewrite of compute_edit_distance.py, whose per-character CIGAR walk and
// repeated per-window scans made it slow on whole chromosomes. Output on stdout
// is identical.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "helpers.hpp"

namespace {

constexpr int64_t indel_length_cutoff = 25;

// Mark the inverted intervals listed in a BED-like file. Absent files are
// simply skipped, matching the Python.
void read_inversions(const std::string &fname,
                     const std::string &ref_header, const std::string &qry_header,
                     std::vector<char> *ref_inv, std::vector<char> *qry_inv) {
    std::ifstream in(fname);
    if (!in)
        return;

    std::string line;
    while (std::getline(in, line)) {
        const std::vector<std::string_view> fields = split(line, '\t');
        if (fields.size() < 9)
            continue;

        if (fields[0] != ref_header)
            throw std::runtime_error(fname + ": reference name does not match the FASTA");
        if (fields[6] != qry_header)
            throw std::runtime_error(fname + ": query name does not match the FASTA");

        const int64_t ref_begin = to_int(fields[1]);
        const int64_t ref_end = to_int(fields[2]);
        const int64_t qry_begin = to_int(fields[7]);
        const int64_t qry_end = to_int(fields[8]);

        std::fill(ref_inv->begin() + ref_begin, ref_inv->begin() + ref_end, char(1));
        std::fill(qry_inv->begin() + qry_begin, qry_inv->begin() + qry_end, char(1));
    }
}

bool same_inversion_state(const std::vector<char> &ref_inv, int64_t ref_pos,
                          const std::vector<char> &qry_inv, int64_t qry_pos,
                          int64_t length) {
    for (int64_t i = 0; i < length; ++i) {
        if (ref_inv[ref_pos + i] != qry_inv[qry_pos + i])
            return false;
    }
    return true;
}

// Everything one aligned window contributes.
struct Window {
    int64_t either_N = 0;   // positions where either sequence has an N
    int64_t observed_eq = 0;
    int64_t observed_neq = 0;
    int64_t ref_N = 0;
    int64_t qry_N = 0;
};

Window scan(std::string_view ref, std::string_view qry) {
    Window w;
    for (size_t i = 0; i < ref.size(); ++i) {
        const char r = ref[i];
        const char q = qry[i];
        const bool r_is_N = (r == 'N');
        const bool q_is_N = (q == 'N');
        w.ref_N += r_is_N;
        w.qry_N += q_is_N;
        if (r_is_N || q_is_N) {
            ++w.either_N;
        } else if (r == q) {
            ++w.observed_eq;
        } else {
            ++w.observed_neq;
        }
    }
    return w;
}

struct Totals {
    int64_t length_a = 0;
    int64_t length_b = 0;
    int64_t full_length_a = 0;
    int64_t full_length_b = 0;
    int64_t dist = 0;
    int64_t id = 0;
    int64_t matches = 0;
    int64_t length_short_indels_a = 0;
    int64_t length_short_indels_b = 0;
    int64_t nchrom = 0;
};

void process(const std::string &ref_fname, const std::string &qry_fname,
             const std::string &cigar_fname, const std::string &inv_fname,
             bool check_N, bool check_eq, Totals *totals) {
    // name is the CIGAR file's basename up to the first dot, minus "chr"
    std::string name = cigar_fname.substr(cigar_fname.find_last_of('/') + 1);
    name = name.substr(0, name.find('.'));
    std::cerr << name << std::endl;
    const std::string label = name.size() > 3 ? name.substr(3) : std::string();

    std::string ref_header, qry_header;
    const std::string ref_seq = read_fasta(ref_fname, &ref_header);
    const std::string qry_seq = read_fasta(qry_fname, &qry_header);

    std::vector<char> ref_inv(ref_seq.size(), 0);
    std::vector<char> qry_inv(qry_seq.size(), 0);
    read_inversions(inv_fname, ref_header, qry_header, &ref_inv, &qry_inv);

    const std::string cigar = read_cigar(cigar_fname);
    if (cigar.empty())
        return;

    int64_t length = 0, length_q = 0;
    int64_t short_indels_a = 0, short_indels_b = 0;
    int64_t dist = 0, cur_id = 0, cur_matches = 0;
    int64_t cur_diag = 0, max_diag = 0;
    int64_t ref_N = 0, qry_N = 0;
    int64_t cur_num = 0;

    // Whether the aligned runs are measured against the sequences or believed.
    const bool verify = check_N || check_eq;

    for (const char c : cigar) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            cur_num = cur_num * 10 + (c - '0');
            continue;
        }

        // X, D, I, S and H all count towards the edit distance
        if (c == 'X' || c == 'D' || c == 'I' || c == 'S' || c == 'H')
            dist += cur_num;

        if (c == 'I') {
            cur_diag -= cur_num;
            const int64_t n = count_N(std::string_view(qry_seq).substr(length_q, cur_num));
            if (std::abs(cur_diag) > std::abs(max_diag))
                max_diag = cur_diag;
            if (cur_num - n <= indel_length_cutoff)
                short_indels_b += cur_num - n;
            qry_N += n;
            length_q += cur_num;
        } else if (c == 'D') {
            cur_diag += cur_num;
            const int64_t n = count_N(std::string_view(ref_seq).substr(length, cur_num));
            if (std::abs(cur_diag) > std::abs(max_diag))
                max_diag = cur_diag;
            if (cur_num - n <= indel_length_cutoff)
                short_indels_a += cur_num - n;
            ref_N += n;
            length += cur_num;
        } else if (c == 'M') {
            // M marks positions aligned against an N, so every base in the run
            // should have an N on at least one side; = and X only ever cover
            // non-N positions. It consumes both sequences like = and X do.
            //
            // Believed, the run still has to be looked at once to split its N
            // count between the two sides, the same way D and I runs are. Under
            // either check the label is not taken on trust: a position with no
            // N on either side belongs in an = or an X run, so it is reported
            // and rejected, and the run is scored by what its bases actually
            // are, leaving N positions out of the identity as before but
            // counting any real match or mismatch hiding inside it.
            if (!same_inversion_state(ref_inv, length, qry_inv, length_q, cur_num))
                throw std::runtime_error("M run spans a mismatched inversion state");
            const std::string_view ref_w = std::string_view(ref_seq).substr(length, cur_num);
            const std::string_view qry_w = std::string_view(qry_seq).substr(length_q, cur_num);

            if (!verify) {
                ref_N += count_N(ref_w);
                qry_N += count_N(qry_w);
            } else {
                const Window w = scan(ref_w, qry_w);
                if (w.either_N != cur_num) {
                    std::cerr << (length + 1) << "-" << (length + cur_num) << "\n"
                              << (length_q + 1) << "-" << (length_q + cur_num) << "\n"
                              << c << "\t" << "\t=" << " " << w.observed_eq
                              << " X " << w.observed_neq << " N " << w.either_N
                              << " total " << cur_num << "\n"
                              << ref_w << "\n" << qry_w << "\n" << std::endl;
                    throw std::runtime_error("M run contains a position without an N");
                }
                dist += w.observed_neq;
                cur_id += w.observed_eq;
                cur_matches += w.observed_eq + w.observed_neq;
                ref_N += w.ref_N;
                qry_N += w.qry_N;
            }
            length += cur_num;
            length_q += cur_num;
        } else if (c == '=') {
            if (!same_inversion_state(ref_inv, length, qry_inv, length_q, cur_num))
                throw std::runtime_error("= run spans a mismatched inversion state");

            if (!verify) {
                // believed to be a match between two non-N bases
                cur_id += cur_num;
                cur_matches += cur_num;
            } else {
                const std::string_view ref_w = std::string_view(ref_seq).substr(length, cur_num);
                const std::string_view qry_w = std::string_view(qry_seq).substr(length_q, cur_num);
                const Window w = scan(ref_w, qry_w);

                if (w.observed_neq > 0) {
                    std::cerr << (length + 1) << "-" << (length + cur_num) << "\n"
                              << (length_q + 1) << "-" << (length_q + cur_num) << "\n"
                              << c << "\t" << "\t=" << " " << w.observed_eq
                              << " X " << w.observed_neq << " N " << w.either_N
                              << " total " << cur_num << "\n"
                              << ref_w << "\n" << qry_w << "\n" << std::endl;
                }
                if (check_eq && w.observed_neq != 0)
                    throw std::runtime_error("= run contains differing bases");

                cur_id += w.observed_eq;
                cur_matches += w.observed_eq;
                ref_N += w.ref_N;
                qry_N += w.qry_N;
            }
            length += cur_num;
            length_q += cur_num;
        } else if (c == 'X') {
            if (!same_inversion_state(ref_inv, length, qry_inv, length_q, cur_num))
                throw std::runtime_error("X run spans a mismatched inversion state");

            if (!verify) {
                // believed to be a mismatch between two non-N bases
                cur_matches += cur_num;
            } else {
                const Window w = scan(std::string_view(ref_seq).substr(length, cur_num),
                                      std::string_view(qry_seq).substr(length_q, cur_num));
                if (check_eq && w.observed_eq != 0)
                    throw std::runtime_error("X run contains identical bases");

                cur_matches += w.observed_neq - w.either_N;
                ref_N += w.ref_N;
                qry_N += w.qry_N;
            }
            length += cur_num;
            length_q += cur_num;
        }

        if (length > int64_t(ref_seq.size()))
            throw std::runtime_error("CIGAR consumes more reference than exists");
        if (length_q > int64_t(qry_seq.size()))
            throw std::runtime_error("CIGAR consumes more query than exists");
        if (cur_diag != length - length_q)
            throw std::runtime_error("diagonal is out of step with the consumed lengths");
        cur_num = 0;
    }

    if (length != int64_t(ref_seq.size()))
        throw std::runtime_error("CIGAR does not span the whole reference");
    if (length_q != int64_t(qry_seq.size()))
        throw std::runtime_error("CIGAR does not span the whole query");

    const int64_t full_length = length;
    const int64_t full_length_q = length_q;
    length -= ref_N;
    length_q -= qry_N;

    const int64_t minlen = std::min(length, length_q);
    const int64_t maxlen = std::max(length, length_q);
    std::printf("%s\t%lld->%lld\t%lld->%lld\t%lld\t%.2f\t%.2f\n",
                label.c_str(),
                (long long)full_length, (long long)length,
                (long long)full_length_q, (long long)length_q,
                (long long)cur_id,
                100.0 * double(cur_id) / double(maxlen),
                100.0 * double(cur_id) / double(minlen));

    totals->full_length_a += full_length;
    totals->full_length_b += full_length_q;
    totals->dist += dist;
    totals->length_a += length;
    totals->length_b += length_q;
    totals->id += cur_id;
    totals->matches += cur_matches;
    totals->length_short_indels_a += cur_matches + short_indels_a;
    totals->length_short_indels_b += cur_matches + short_indels_b;
    ++totals->nchrom;
}

}  // namespace

int main(int argc, char **argv) {
    if (argc < 3 || (argc - 3) % 4 != 0) {
        std::cerr << "usage: " << argv[0]
                  << " CHECK_N CHECK_EQ (REF.fa QRY.fa CIGAR INV.bed)...\n";
        return 2;
    }

    const bool check_N = std::atoi(argv[1]) != 0;
    const bool check_eq = std::atoi(argv[2]) != 0;
    std::cerr << (check_eq ? "True" : "False") << " "
              << (check_N ? "True" : "False") << std::endl;

    try {
        Totals totals;
        for (int i = 3; i + 3 < argc; i += 4)
            process(argv[i], argv[i+1], argv[i+2], argv[i+3], check_N, check_eq, &totals);

        if (totals.nchrom > 1) {
            const int64_t minlen = std::min(totals.length_a, totals.length_b);
            const int64_t maxlen = std::max(totals.length_a, totals.length_b);
            std::printf("HG002\t%lld->%lld\t%lld->%lld\t%lld\t%.2f\t%.2f\n",
                        (long long)totals.full_length_a, (long long)totals.length_a,
                        (long long)totals.full_length_b, (long long)totals.length_b,
                        (long long)totals.id,
                        100.0 * double(totals.id) / double(maxlen),
                        100.0 * double(totals.id) / double(minlen));
        }
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
