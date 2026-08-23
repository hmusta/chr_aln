// Summarise a global CIGAR alignment against the two sequences it aligns.
//
// Usage:
//     compute_edit_distance CHECK_EQ (REF.fa QRY.fa CIGAR INV.bed)...
//
// CHECK_EQ is 0/1. Each group of four arguments describes one chromosome;
// INV.bed may be absent, in which case it is skipped.
//
// The CIGAR uses = and X only for positions where neither sequence has an N,
// and M for positions aligned against an N on either side. All three consume
// reference and query alike; D and I consume one side each.
//
// With CHECK_EQ clear the aligned operations are taken at their word: = is a
// match between two non-N bases, X a mismatch between two non-N bases, and M a
// position against an N. Nothing is read back from the sequences except the N
// counts inside M, D and I runs, which are needed to shorten the lengths.
//
// With CHECK_EQ set nothing is assumed. Every =, X and M run is re-read and
// each position counted as whatever it actually is, so a run carrying the wrong
// operation is corrected rather than believed. This is what a CIGAR with no M
// convention needs, since there the Ns hide inside = and X runs and can only be
// found by looking. The number of corrected positions is reported on stderr.
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
             bool check_eq, Totals *totals) {
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
    int64_t corrected = 0;

    for (const char c : cigar) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            cur_num = cur_num * 10 + (c - '0');
            continue;
        }

        // D, I, S and H always count towards the edit distance; what an
        // aligned run contributes depends on whether it is believed or re-read
        if (c == 'D' || c == 'I' || c == 'S' || c == 'H')
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
        } else if (c == 'M' || c == '=' || c == 'X') {
            if (!same_inversion_state(ref_inv, length, qry_inv, length_q, cur_num))
                throw std::runtime_error(std::string(1, c)
                                         + " run spans a mismatched inversion state");
            const std::string_view ref_w = std::string_view(ref_seq).substr(length, cur_num);
            const std::string_view qry_w = std::string_view(qry_seq).substr(length_q, cur_num);

            if (!check_eq) {
                // Believed. Only M is looked at, and only to split its N count
                // between the two sides so the lengths can be shortened; = and
                // X are non-N by convention and contribute no N at all.
                if (c == '=') {
                    cur_id += cur_num;
                    cur_matches += cur_num;
                } else if (c == 'X') {
                    dist += cur_num;
                    cur_matches += cur_num;
                } else {
                    ref_N += count_N(ref_w);
                    qry_N += count_N(qry_w);
                }
            } else {
                // Re-read. Each position counts as whatever it actually is: an
                // N on either side keeps it out of the identity, two equal
                // bases make it a match, two differing bases a mismatch. A run
                // carrying the wrong operation is corrected, not rejected.
                const Window w = scan(ref_w, qry_w);

                if (c == '=' && w.observed_neq > 0) {
                    std::cerr << (length + 1) << "-" << (length + cur_num) << "\n"
                              << (length_q + 1) << "-" << (length_q + cur_num) << "\n"
                              << c << "\t" << "\t=" << " " << w.observed_eq
                              << " X " << w.observed_neq << " N " << w.either_N
                              << " total " << cur_num << "\n"
                              << ref_w << "\n" << qry_w << "\n" << std::endl;
                }

                // positions that really were what the operation claimed
                const int64_t as_labelled = (c == '=') ? w.observed_eq
                                          : (c == 'X') ? w.observed_neq
                                                       : w.either_N;
                corrected += cur_num - as_labelled;

                dist += w.observed_neq;
                cur_id += w.observed_eq;
                cur_matches += w.observed_eq + w.observed_neq;
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

    if (corrected > 0)
        std::cerr << name << ": corrected " << corrected
                  << " position(s) whose operation disagreed with the sequences"
                  << std::endl;

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
    if (argc < 2 || (argc - 2) % 4 != 0) {
        std::cerr << "usage: " << argv[0]
                  << " CHECK_EQ (REF.fa QRY.fa CIGAR INV.bed)...\n";
        return 2;
    }

    const bool check_eq = std::atoi(argv[1]) != 0;
    std::cerr << (check_eq ? "True" : "False") << std::endl;

    try {
        Totals totals;
        for (int i = 2; i + 3 < argc; i += 4)
            process(argv[i], argv[i+1], argv[i+2], argv[i+3], check_eq, &totals);

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
