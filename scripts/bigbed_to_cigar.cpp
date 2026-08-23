// Build a global CIGAR alignment for each aligned chromosome pair from a
// bigBed of heterozygous sites.
//
// Usage:
//     bigbed_to_cigar [options] HETSITES.bb OUT_PREFIX
//
// The bigBed lists every position at which the two haplotypes differ. Each
// record sits on one haplotype and names its mate:
//
//     chr22_PATERNAL  3036 3037  chr22_MATERNAL_2167_C_T_F
//
// so a record pairs a reference interval with a query interval. Between two
// consecutive records the haplotypes are identical, which is what makes a
// CIGAR derivable: the stretch between sites becomes '=' and the site itself
// becomes 'X', 'I' or 'D'. Where the two sides advance by different amounts
// the file says nothing about how they align, and that region is emitted as
// D/I rather than guessed at.
//
// Operations, matching the other tools here:
//     =   identical bases, neither an N, consumes reference and query
//     X   differing bases, neither an N, consumes reference and query
//     M   a position where either side is an N, consumes reference and query
//     D   consumes reference only
//     I   consumes query only
//
// Sites inside an inversion are flagged 'R' and run backwards along the query.
// A CIGAR is colinear, so those regions are reverse-complemented in the query
// first, which makes them ascend again and lets them carry matches. The
// rewritten query is written next to the CIGAR, since the CIGAR describes the
// alignment against that sequence rather than the original. Detecting the
// inversions needs no FASTA, but rewriting the query does: without --fasta-dir
// the inverted runs are left unaligned as D/I instead.

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <unordered_map>

#include <cstring>

#include <zlib.h>

#include "helpers.hpp"

namespace {

constexpr uint32_t BIGBED_MAGIC = 0x8789F2EB;
constexpr uint32_t BPT_MAGIC    = 0x78CA8C91;
constexpr uint32_t RTREE_MAGIC  = 0x2468ACE0;

// ------------------------------------------------------------ raw reading --

struct File {
    std::ifstream in;

    explicit File(const std::string &fname) : in(fname, std::ios::binary) {
        if (!in)
            throw std::runtime_error("cannot open " + fname);
    }

    std::string read(uint64_t offset, size_t length) {
        in.seekg(static_cast<std::streamoff>(offset));
        std::string buf(length, '\0');
        if (!in.read(buf.data(), static_cast<std::streamsize>(length)))
            throw std::runtime_error("short read in the bigBed");
        return buf;
    }
};

// bigBed is written in the host order of whichever machine created it; the
// magic tells us which, and every file in practice is little-endian.
template <class T>
T get(const std::string &buf, size_t offset) {
    T value;
    if (offset + sizeof(T) > buf.size())
        throw std::runtime_error("truncated bigBed structure");
    std::memcpy(&value, buf.data() + offset, sizeof(T));
    return value;
}

std::string inflate_block(const std::string &raw, size_t hint) {
    std::string out(std::max<size_t>(hint, raw.size() * 4), '\0');
    while (true) {
        uLongf size = static_cast<uLongf>(out.size());
        const int rc = uncompress(reinterpret_cast<Bytef *>(out.data()), &size,
                                  reinterpret_cast<const Bytef *>(raw.data()),
                                  static_cast<uLong>(raw.size()));
        if (rc == Z_OK) {
            out.resize(size);
            return out;
        }
        if (rc != Z_BUF_ERROR)
            throw std::runtime_error("zlib failed to expand a data block");
        out.resize(out.size() * 2);
    }
}

// ---------------------------------------------------------------- bigBed --

struct Chrom {
    std::string name;
    uint32_t id = 0;
    uint32_t size = 0;
};

struct BigBed {
    File file;
    uint64_t chrom_tree_offset = 0;
    uint64_t full_index_offset = 0;
    uint32_t uncompress_buf_size = 0;
    std::vector<Chrom> chroms;                  // by id
    std::unordered_map<std::string, uint32_t> by_name;

    explicit BigBed(const std::string &fname) : file(fname) {
        const std::string hdr = file.read(0, 64);
        if (get<uint32_t>(hdr, 0) != BIGBED_MAGIC)
            throw std::runtime_error(fname + ": not a bigBed file");

        chrom_tree_offset = get<uint64_t>(hdr, 8);
        full_index_offset = get<uint64_t>(hdr, 24);
        uncompress_buf_size = get<uint32_t>(hdr, 52);

        read_chrom_tree();
    }

    void read_chrom_tree() {
        const std::string hdr = file.read(chrom_tree_offset, 32);
        if (get<uint32_t>(hdr, 0) != BPT_MAGIC)
            throw std::runtime_error("bad chromosome B+ tree");
        const uint32_t key_size = get<uint32_t>(hdr, 8);
        const uint32_t val_size = get<uint32_t>(hdr, 12);
        if (val_size != 8)
            throw std::runtime_error("unexpected chromosome record size");

        std::vector<Chrom> found;
        walk_bpt(chrom_tree_offset + 32, key_size, val_size, found);

        uint32_t max_id = 0;
        for (const Chrom &c : found)
            max_id = std::max(max_id, c.id);
        chroms.assign(max_id + 1, Chrom{});
        for (const Chrom &c : found) {
            chroms[c.id] = c;
            by_name[c.name] = c.id;
        }
    }

    void walk_bpt(uint64_t offset, uint32_t key_size, uint32_t val_size,
                  std::vector<Chrom> &out) {
        const std::string node = file.read(offset, 4);
        const bool is_leaf = get<uint8_t>(node, 0) != 0;
        const uint16_t count = get<uint16_t>(node, 2);

        const size_t item = key_size + (is_leaf ? val_size : 8);
        const std::string body = file.read(offset + 4, count * item);

        if (is_leaf) {
            for (uint16_t i = 0; i < count; ++i) {
                const size_t at = i * item;
                std::string key = body.substr(at, key_size);
                key.resize(std::min(key.size(), key.find('\0') == std::string::npos
                                                    ? key.size() : key.find('\0')));
                Chrom c;
                c.name = key;
                c.id = get<uint32_t>(body, at + key_size);
                c.size = get<uint32_t>(body, at + key_size + 4);
                out.push_back(std::move(c));
            }
        } else {
            for (uint16_t i = 0; i < count; ++i)
                walk_bpt(get<uint64_t>(body, i * item + key_size), key_size, val_size, out);
        }
    }

    struct Block {
        uint32_t start_chrom, end_chrom;
        uint64_t offset, size;
    };

    std::vector<Block> leaf_blocks() {
        const std::string hdr = file.read(full_index_offset, 48);
        if (get<uint32_t>(hdr, 0) != RTREE_MAGIC)
            throw std::runtime_error("bad R-tree index");
        std::vector<Block> out;
        walk_rtree(full_index_offset + 48, out);
        return out;
    }

    void walk_rtree(uint64_t offset, std::vector<Block> &out) {
        const std::string node = file.read(offset, 4);
        const bool is_leaf = get<uint8_t>(node, 0) != 0;
        const uint16_t count = get<uint16_t>(node, 2);
        const size_t item = is_leaf ? 32 : 24;
        const std::string body = file.read(offset + 4, count * item);

        for (uint16_t i = 0; i < count; ++i) {
            const size_t at = i * item;
            if (is_leaf) {
                out.push_back({ get<uint32_t>(body, at), get<uint32_t>(body, at + 8),
                                get<uint64_t>(body, at + 16), get<uint64_t>(body, at + 24) });
            } else {
                walk_rtree(get<uint64_t>(body, at + 16), out);
            }
        }
    }

    // Visit every BED record, decompressing blocks as needed.
    template <class Visit>
    void for_each_record(Visit &&visit) {
        for (const Block &b : leaf_blocks()) {
            const std::string raw = file.read(b.offset, b.size);
            const std::string buf = uncompress_buf_size
                                        ? inflate_block(raw, uncompress_buf_size)
                                        : raw;
            size_t p = 0;
            while (p + 12 <= buf.size()) {
                const uint32_t chrom_id = get<uint32_t>(buf, p);
                const uint32_t start = get<uint32_t>(buf, p + 4);
                const uint32_t end = get<uint32_t>(buf, p + 8);
                p += 12;
                const size_t nul = buf.find('\0', p);
                if (nul == std::string::npos)
                    throw std::runtime_error("unterminated name in a data block");
                visit(chrom_id, start, end, std::string_view(buf.data() + p, nul - p));
                p = nul + 1;
            }
        }
    }
};

// ------------------------------------------------------------- het sites --

struct Interval {
    int64_t begin, end;
};

struct Site {
    int64_t rbegin = 0;     // reference interval, from the record's own coords
    int64_t rend = 0;
    int64_t qbegin = 0;     // query start, from the name (1-based there)
    int64_t qlen = -1;      // query allele length; -1 until recovered
    bool inverted = false;
    bool complete = false;  // the name carried both alleles and the strand
};

// name is <mate>_<pos>_<thisAllele>_<otherAllele>_<F|R>, and is truncated at
// 255 characters, so a long allele can lose its tail. The position always
// survives; the query length is recovered from the mate record when it does not.
bool parse_name(std::string_view name, std::string_view mate, Site *site) {
    if (name.size() < mate.size() + 2 || name.compare(0, mate.size(), mate) != 0
        || name[mate.size()] != '_')
        return false;

    size_t p = mate.size() + 1;
    const size_t digits = p;
    while (p < name.size() && std::isdigit(static_cast<unsigned char>(name[p])))
        ++p;
    if (p == digits || p >= name.size() || name[p] != '_')
        return false;
    site->qbegin = to_int(name.substr(digits, p - digits)) - 1;   // 1-based in the name
    ++p;

    const size_t a1 = name.find('_', p);
    if (a1 == std::string_view::npos)
        return true;                       // truncated inside the first allele
    const size_t a2 = name.find('_', a1 + 1);
    if (a2 == std::string_view::npos)
        return true;                       // truncated inside the second allele

    const std::string_view other = name.substr(a1 + 1, a2 - a1 - 1);
    const std::string_view strand = name.substr(a2 + 1);
    if (strand != "F" && strand != "R")
        return true;

    site->qlen = (other == "*") ? 0 : static_cast<int64_t>(other.size());
    site->inverted = (strand == "R");
    site->complete = true;
    return true;
}

struct PairStats {
    int64_t sites = 0, dropped = 0, recovered = 0, unrecovered = 0;
    int64_t drop_unresolved = 0, drop_ref_back = 0, drop_qry_back = 0;
    int64_t inverted_sites = 0, inverted_runs = 0;
    int64_t flipped_regions = 0, flipped_bp = 0, flip_conflicts = 0;
    int64_t unequal_gaps = 0, unaligned_ref = 0, unaligned_qry = 0;
};

// Group the inverted sites into inverted regions.
//
// A single inversion is not a single run: forward-flagged sites are scattered
// through it and would split it into fragments, and reverse-complementing each
// fragment in place is not the same operation as flipping the whole region.
// What actually marks one inversion is that the query start keeps descending
// as the reference advances, so the sites are clustered on that instead, and a
// cluster ends where the query jumps forward again.
std::vector<Interval> inverted_regions(const std::vector<Site> &sorted_sites) {
    std::vector<Interval> regions;
    bool open = false;
    int64_t prev_q = 0;
    Interval cur{0, 0};

    for (const Site &s : sorted_sites) {
        if (!s.inverted || s.qlen < 0)
            continue;
        const int64_t begin = s.qbegin;
        const int64_t end = s.qbegin + s.qlen;

        if (open && begin < prev_q) {
            cur.begin = std::min(cur.begin, begin);
            cur.end = std::max(cur.end, end);
        } else {
            if (open)
                regions.push_back(cur);
            cur = {begin, end};
            open = true;
        }
        prev_q = begin;
    }
    if (open)
        regions.push_back(cur);

    // Overlapping regions cannot each be reverse-complemented, so merge them.
    std::sort(regions.begin(), regions.end(),
              [](const Interval &a, const Interval &b) { return a.begin < b.begin; });
    std::vector<Interval> merged;
    for (const Interval &r : regions) {
        if (!merged.empty() && r.begin < merged.back().end)
            merged.back().end = std::max(merged.back().end, r.end);
        else
            merged.push_back(r);
    }
    return merged;
}

// Reverse-complement each region of the query, rebuilding the sequence.
std::string invert_query(const std::string &qry, const std::vector<Interval> &regions) {
    std::string out;
    out.reserve(qry.size());
    int64_t prev = 0;
    for (const Interval &r : regions) {
        out.append(qry, prev, r.begin - prev);
        for (int64_t i = r.end - 1; i >= r.begin; --i)
            out.push_back(complement(qry[i]));
        prev = r.end;
    }
    out.append(qry, prev, qry.size() - prev);
    return out;
}

// Map a query interval through the reverse-complementing of its region.
bool remap_query(int64_t *begin, int64_t *end, const std::vector<Interval> &regions) {
    for (const Interval &r : regions) {
        if (*begin >= r.begin && *end <= r.end) {
            const int64_t b = r.begin + r.end - *end;
            const int64_t e = r.begin + r.end - *begin;
            *begin = b;
            *end = e;
            return true;
        }
        // an interval straddling a boundary has no consistent image
        if (*begin < r.end && *end > r.begin)
            return false;
    }
    return true;
}

// Walk the sites in reference order, emitting '=' between them and the
// difference itself at each one.
std::vector<Op> build_cigar(std::vector<Site> &sites, int64_t ref_len, int64_t qry_len,
                            PairStats *stats) {
    std::vector<Op> cigar;
    int64_t r = 0, q = 0;
    bool in_inverted_run = false;

    for (const Site &s : sites) {
        if (s.qlen < 0 || s.rbegin < r || s.qbegin < q) {
            ++stats->dropped;             // unrecoverable or out of order
            if (s.qlen < 0)          ++stats->drop_unresolved;
            else if (s.rbegin < r)   ++stats->drop_ref_back;
            else                     ++stats->drop_qry_back;
            continue;
        }

        if (s.inverted) {
            // A colinear CIGAR cannot hold an inversion, so the run is left
            // unaligned; the gap accounting below absorbs it.
            if (!in_inverted_run) {
                ++stats->inverted_runs;
                in_inverted_run = true;
            }
            ++stats->inverted_sites;
            continue;
        }
        in_inverted_run = false;

        const int64_t dr = s.rbegin - r;
        const int64_t dq = s.qbegin - q;
        if (dr == dq) {
            append_op(cigar, '=', dr);
        } else {
            // the file does not say how these align, so do not invent one
            ++stats->unequal_gaps;
            stats->unaligned_ref += dr;
            stats->unaligned_qry += dq;
            append_op(cigar, 'D', dr);
            append_op(cigar, 'I', dq);
        }
        r = s.rbegin;
        q = s.qbegin;

        const int64_t rl = s.rend - s.rbegin;
        append_op(cigar, 'X', std::min(rl, s.qlen));
        append_op(cigar, 'D', rl - s.qlen);
        append_op(cigar, 'I', s.qlen - rl);
        r += rl;
        q += s.qlen;
        ++stats->sites;
    }

    const int64_t dr = ref_len - r;
    const int64_t dq = qry_len - q;
    if (dr == dq) {
        append_op(cigar, '=', dr);
    } else {
        if (dr || dq) {
            ++stats->unequal_gaps;
            stats->unaligned_ref += dr;
            stats->unaligned_qry += dq;
        }
        append_op(cigar, 'D', dr);
        append_op(cigar, 'I', dq);
    }
    return cigar;
}

// Relabel the operations that consume both sequences against the sequences
// themselves, so '=' and 'X' never cover an N.
std::vector<Op> label_operations(const std::vector<Op> &cigar,
                                 const std::string &ref, const std::string &qry) {
    std::vector<Op> out;
    int64_t r = 0, q = 0;
    for (const Op &o : cigar) {
        if (!aligns_both(o.op)) {
            append_op(out, o.op, o.length);
            if (consumes_ref(o.op)) r += o.length;
            if (consumes_qry(o.op)) q += o.length;
            continue;
        }
        for (int64_t i = 0; i < o.length; ++i) {
            const char a = ref[r + i], b = qry[q + i];
            append_op(out, (a == 'N' || b == 'N') ? 'M' : (a == b ? '=' : 'X'), 1);
        }
        r += o.length;
        q += o.length;
    }
    return out;
}

struct Verification {
    int64_t ref_consumed = 0, qry_consumed = 0;
    int64_t bad_eq = 0, bad_x = 0, bad_m = 0;
};

Verification verify(const std::vector<Op> &cigar, const std::string &ref,
                    const std::string &qry) {
    Verification v;
    int64_t r = 0, q = 0;
    for (const Op &o : cigar) {
        if (aligns_both(o.op)) {
            for (int64_t i = 0; i < o.length; ++i) {
                const char a = ref[r + i], b = qry[q + i];
                const bool has_n = (a == 'N' || b == 'N');
                if (o.op == '=' && (has_n || a != b)) ++v.bad_eq;
                if (o.op == 'X' && (has_n || a == b)) ++v.bad_x;
                if (o.op == 'M' && !has_n)            ++v.bad_m;
            }
            r += o.length;
            q += o.length;
        } else if (o.op == 'D') {
            r += o.length;
        } else if (o.op == 'I') {
            q += o.length;
        }
    }
    v.ref_consumed = r;
    v.qry_consumed = q;
    return v;
}

const char *USAGE =
    "usage: bigbed_to_cigar [options] HETSITES.bb OUT_PREFIX\n"
    "\n"
    "Build one global CIGAR per aligned chromosome pair from a bigBed of\n"
    "heterozygous sites between two haplotypes.\n"
    "\n"
    "positional arguments:\n"
    "  HETSITES.bb   bigBed whose name field is\n"
    "                <mateChrom>_<pos>_<thisAllele>_<otherAllele>_<F|R>\n"
    "  OUT_PREFIX    prefix for the output files, including any directory\n"
    "\n"
    "options:\n"
    "  --ref-suffix S   haplotype used as the reference (default PATERNAL)\n"
    "  --qry-suffix S   haplotype used as the query (default MATERNAL)\n"
    "  --chrom NAME     only this reference chromosome; repeatable\n"
    "  --fasta-dir DIR  verify against DIR/<chrom>.fa and label N positions as M\n"
    "  -h, --help       show this message\n"
    "\n"
    "output:\n"
    "  OUT_PREFIX.<refChrom>.cigar   one CIGAR per pair\n"
    "  A per-pair summary is written to standard error.\n"
    "\n"
    "operations:\n"
    "  =  identical bases, neither an N   X  differing bases, neither an N\n"
    "  M  a position where either side is an N (only with --fasta-dir)\n"
    "  D  consumes reference only         I  consumes query only\n";

}  // namespace

int main(int argc, char **argv) {
    std::vector<std::string> args;
    std::string ref_suffix = "PATERNAL", qry_suffix = "MATERNAL", fasta_dir;
    std::vector<std::string> only;
    bool want_help = false;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help") want_help = true;
        else if (a == "--ref-suffix" && i + 1 < argc) ref_suffix = argv[++i];
        else if (a == "--qry-suffix" && i + 1 < argc) qry_suffix = argv[++i];
        else if (a == "--fasta-dir" && i + 1 < argc) fasta_dir = argv[++i];
        else if (a == "--chrom" && i + 1 < argc) only.push_back(argv[++i]);
        else args.push_back(a);
    }
    if (want_help || args.size() != 2) {
        std::fputs(USAGE, stderr);
        return want_help ? 0 : 2;
    }

    try {
        BigBed bb(args[0]);

        // pair each reference chromosome with its mate
        std::vector<std::pair<uint32_t, uint32_t>> pairs;
        for (const Chrom &c : bb.chroms) {
            if (c.name.size() <= ref_suffix.size()
                || c.name.compare(c.name.size() - ref_suffix.size(), ref_suffix.size(),
                                  ref_suffix) != 0)
                continue;
            if (!only.empty()
                && std::find(only.begin(), only.end(), c.name) == only.end())
                continue;
            const std::string mate =
                c.name.substr(0, c.name.size() - ref_suffix.size()) + qry_suffix;
            auto it = bb.by_name.find(mate);
            if (it == bb.by_name.end()) {
                std::cerr << "WARNING: no mate " << mate << " for " << c.name << "\n";
                continue;
            }
            pairs.emplace_back(c.id, it->second);
        }
        if (pairs.empty())
            throw std::runtime_error("no chromosome pairs matched the suffixes");

        // one pass over the file collects every record we need
        std::unordered_map<uint32_t, std::vector<Site>> ref_sites;
        struct Mate { int64_t length; bool inverted; };
        std::unordered_map<uint32_t, std::unordered_map<int64_t, Mate>> mate_len;
        std::unordered_map<uint32_t, uint32_t> is_ref, is_qry;
        for (const auto &[rid, qid] : pairs) {
            is_ref[rid] = qid;
            is_qry[qid] = rid;
        }

        bb.for_each_record([&](uint32_t chrom_id, uint32_t start, uint32_t end,
                               std::string_view name) {
            if (auto it = is_ref.find(chrom_id); it != is_ref.end()) {
                Site s;
                s.rbegin = start;
                s.rend = end;
                if (parse_name(name, bb.chroms[it->second].name, &s))
                    ref_sites[chrom_id].push_back(s);
            }
            if (auto it = is_qry.find(chrom_id); it != is_qry.end()) {
                Site m;
                parse_name(name, bb.chroms[it->second].name, &m);
                mate_len[chrom_id][start] = { static_cast<int64_t>(end) - start, m.inverted };
            }
        });

        for (const auto &[rid, qid] : pairs) {
            const Chrom &rc = bb.chroms[rid];
            const Chrom &qc = bb.chroms[qid];
            std::vector<Site> &sites = ref_sites[rid];
            const auto &lens = mate_len[qid];

            PairStats stats;
            for (Site &s : sites) {
                if (s.qlen >= 0)
                    continue;
                if (auto it = lens.find(s.qbegin); it != lens.end()) {
                    s.qlen = it->second.length;
                    s.inverted = it->second.inverted;
                    ++stats.recovered;
                } else {
                    ++stats.unrecovered;
                }
            }

            std::sort(sites.begin(), sites.end(), [](const Site &a, const Site &b) {
                return std::tie(a.rbegin, a.rend) < std::tie(b.rbegin, b.rend);
            });

            std::string ref_seq, qry_seq, qry_header;
            std::vector<Interval> regions;
            if (!fasta_dir.empty()) {
                std::string header;
                ref_seq = read_fasta(fasta_dir + "/" + rc.name + ".fa", &header);
                qry_seq = read_fasta(fasta_dir + "/" + qc.name + ".fa", &qry_header);
                if (int64_t(ref_seq.size()) != int64_t(rc.size)
                    || int64_t(qry_seq.size()) != int64_t(qc.size))
                    throw std::runtime_error(rc.name + ": FASTA length disagrees with the bigBed");

                // Reverse-complement the inverted regions so they ascend with the
                // reference, then move every site inside one to its new coordinates.
                regions = inverted_regions(sites);
                if (!regions.empty()) {
                    qry_seq = invert_query(qry_seq, regions);
                    for (Site &s : sites) {
                        if (s.qlen < 0)
                            continue;
                        int64_t b = s.qbegin, e = s.qbegin + s.qlen;
                        if (!remap_query(&b, &e, regions)) {
                            ++stats.flip_conflicts;
                            s.qlen = -1;          // straddles a boundary; drop it
                            continue;
                        }
                        s.qbegin = b;
                        s.inverted = false;       // colinear now
                    }
                    stats.flipped_regions = int64_t(regions.size());
                    for (const Interval &r : regions)
                        stats.flipped_bp += r.end - r.begin;
                }
            }

            std::vector<Op> cigar = build_cigar(sites, rc.size, qc.size, &stats);
            if (!fasta_dir.empty())
                cigar = label_operations(cigar, ref_seq, qry_seq);

            std::map<char, int64_t> counts;
            for (const Op &o : cigar)
                counts[o.op] += o.length;

            std::cerr << rc.name << " vs " << qc.name << "\n"
                      << "  sites used: " << stats.sites
                      << " (recovered from mate: " << stats.recovered
                      << ", unrecoverable: " << stats.unrecovered
                      << ", dropped: " << stats.dropped
                      << " [unresolved " << stats.drop_unresolved
                      << ", ref back " << stats.drop_ref_back
                      << ", qry back " << stats.drop_qry_back << "])\n"
                      << "  inverted regions flipped: " << stats.flipped_regions
                      << " (" << stats.flipped_bp << " bp)"
                      << ", left unaligned: " << stats.inverted_sites
                      << " in " << stats.inverted_runs << " runs"
                      << ", straddling a boundary: " << stats.flip_conflicts << "\n"
                      << "  unequal gaps: " << stats.unequal_gaps
                      << " (ref " << stats.unaligned_ref
                      << ", qry " << stats.unaligned_qry << " unaligned)\n"
                      << "  operations: " << cigar.size();
            for (const char op : std::string("=XMID"))
                if (counts[op]) std::cerr << "  " << op << ": " << counts[op];
            std::cerr << "\n";

            const auto [span_r, span_q] = ops_span(cigar);
            std::cerr << "  reference consumed: " << span_r << " / " << rc.size
                      << (span_r == int64_t(rc.size) ? " OK" : " MISMATCH")
                      << "\n  query consumed:     " << span_q << " / " << qc.size
                      << (span_q == int64_t(qc.size) ? " OK" : " MISMATCH") << "\n";

            if (!fasta_dir.empty()) {
                const Verification v = verify(cigar, ref_seq, qry_seq);
                std::cerr << "  '=' not identical: " << v.bad_eq
                          << "   'X' identical: " << v.bad_x
                          << "   'M' without an N: " << v.bad_m << "\n";
            }
            std::fprintf(stderr, "  identity: %lld / %lld (%.2f%%)\n",
                         (long long)counts['='], (long long)std::min(rc.size, qc.size),
                         100.0 * double(counts['=']) / double(std::min(rc.size, qc.size)));

            const std::string out_name = args[1] + "." + rc.name + ".cigar";
            std::ofstream out(out_name);
            if (!out)
                throw std::runtime_error("cannot write " + out_name);
            for (const Op &o : cigar)
                out << o.length << o.op;
            out << "\n";
            std::cerr << "  wrote " << out_name << "\n";

            if (!fasta_dir.empty()) {
                // the CIGAR describes the alignment against this sequence, which
                // differs from the original query wherever an inversion was flipped
                const std::string qry_name = args[1] + "." + qc.name + ".qry.fa";
                std::ofstream qout(qry_name);
                if (!qout)
                    throw std::runtime_error("cannot write " + qry_name);
                qout << ">" << (regions.empty() ? "" : "inverted_") << qry_header << "\n";
                for (size_t i = 0; i < qry_seq.size(); i += 60)
                    qout << qry_seq.substr(i, 60) << "\n";
                std::cerr << "  wrote " << qry_name << "\n";
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
