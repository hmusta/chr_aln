// Build a global CIGAR alignment between a reference and a query from SyRI output.
//
// Usage:
//     syri_to_cigar [--format FMT] REF.fa QRY.fa SYRI.out ALIGNMENTS OUT_PREFIX
//
// SyRI classifies each alignment (syntenic, inverted, duplicated, translocated),
// which is what decides whether a block can join a colinear chain, but it only
// summarises the alignment itself as variant records. The alignments it was run
// on carry the real base-level alignment, so the classification is taken from
// SyRI and the alignment from minimap2 PAF or nucmer delta. Blocks with no
// matching alignment record fall back to reconstructing one from the variant
// records, which is also what show-coords input leaves every block doing, since
// coords carry no base-level alignment at all.
//
// Operations:
//     =   identical bases, neither an N, consumes reference and query
//     X   differing bases, neither an N, consumes reference and query
//     M   a position where either side is an N, consumes reference and query
//     D   consumes reference only
//     I   consumes query only
//
// The alignment is global: the operations consuming reference sum to the
// reference length and those consuming query sum to the query length.
//
// Inverted regions are reverse-complemented in the query before chaining, which
// makes them colinear with the reference so they contribute matches rather than
// indels. Every query coordinate inside such a region is remapped to match. The
// reverse-complemented query is written next to the CIGAR, since the CIGAR
// describes the alignment against that sequence and not the original query.
//
// A CIGAR is still colinear, so translocations and duplicate copies cannot be
// reached and come out as D/I, as does anything SyRI left unaligned.

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>

#include "helpers.hpp"

namespace {

// blocks that can carry matches once inversions have been flipped
const std::set<std::string> CHAIN_TYPES = {"SYNAL", "INVAL"};
// Blocks whose query span is reverse-complemented. These are the parents of the
// inverted alignments in CHAIN_TYPES, not the alignments themselves: SyRI splits
// one inversion into several alignments whose query intervals can overlap each
// other, and overlapping intervals cannot each be flipped. Flipping the whole
// block instead is equivalent, because reverse-complementing a region also
// reverse-complements every interval inside it, just relocated.
//
// Only true inversions qualify. Inverted translocations and inverted
// duplications are not reachable in a colinear chain, so flipping them would
// clip neighbouring forward blocks for no gain.
const std::set<std::string> INVERTED_TYPES = {"INV"};

enum class VarType { SNP, INS, DEL, HDR };

constexpr int64_t MAX_SHIFT = 4096;

struct Interval {
    int64_t begin, end;
};

using Key = std::array<int64_t, 4>;

struct Alignment {
    char strand = '+';
    bool has_ops = false;
    std::vector<Op> ops;
};

struct Variant {
    int64_t rbegin, rend, qbegin, qend;
    VarType type;
    int parent;
};

struct Block {
    int64_t rbegin, rend, qbegin, qend;
    std::string id;
    std::string type;
    int parent;
    Key key;
};

// ---------------------------------------------------------------- SyRI ------

// SyRI writes one-based inclusive coordinates; these come back zero-based and
// half-open. A '-' start or end marks a record unaligned on that side and
// becomes zero, which is how NOTAL records are recognised.
template <class Visit>
void stream_syri(const std::string &fname, Visit &&visit) {
    std::ifstream in(fname);
    if (!in)
        throw std::runtime_error("cannot open " + fname);

    std::string line;
    int64_t number = 0;
    while (std::getline(in, line)) {
        ++number;
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();
        const std::vector<std::string_view> f = split_ws(line);
        if (f.size() < 11)
            continue;
        try {
            const int64_t rbegin = f[1] == "-" ? 0 : to_int(f[1]) - 1;
            const int64_t rend   = f[2] == "-" ? 0 : to_int(f[2]);
            const int64_t qbegin = f[6] == "-" ? 0 : to_int(f[6]) - 1;
            const int64_t qend   = f[7] == "-" ? 0 : to_int(f[7]);
            if (rbegin < 0 || qbegin < 0)
                throw std::runtime_error("negative coordinate");
            visit(rbegin, rend, qbegin, qend, f[8], f[9], f[10]);
        } catch (const std::exception &e) {
            throw std::runtime_error(fname + ":" + std::to_string(number) + ": " + e.what());
        }
    }
}

// ----------------------------------------------------------- alignments -----

// The PAF is the alignment SyRI was run on, so its records carry the exact
// base-level alignment that SyRI only summarises as variant records. Query
// coordinates are always on the forward strand; for a '-' record the CIGAR
// describes the reverse complement of that interval, which is precisely the
// sequence we get after flipping an inverted region.
std::map<Key, Alignment> read_paf(const std::string &fname, const std::set<Key> &keys) {
    std::ifstream in(fname);
    if (!in)
        throw std::runtime_error("cannot open " + fname);

    std::map<Key, Alignment> alignments;
    std::string line;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();
        const std::vector<std::string_view> f = split(line, '\t');
        if (f.size() < 12)
            continue;

        const Key key = {to_int(f[7]), to_int(f[8]), to_int(f[2]), to_int(f[3])};
        if (!keys.empty() && keys.find(key) == keys.end())
            continue;

        std::string_view cigar;
        bool found = false;
        for (size_t i = 12; i < f.size(); ++i) {
            if (f[i].rfind("cg:Z:", 0) == 0) {
                cigar = f[i].substr(5);
                found = true;
                break;
            }
        }
        if (!found)
            continue;

        Alignment alignment;
        alignment.strand = f[4].empty() ? '+' : f[4][0];
        alignment.has_ops = true;
        int64_t number = 0;
        for (const char c : cigar) {
            if (std::isdigit(static_cast<unsigned char>(c))) {
                number = number * 10 + (c - '0');
            } else {
                append_op(alignment.ops, c, number);
                number = 0;
            }
        }
        alignments[key] = std::move(alignment);
    }
    return alignments;
}

// Each delta value counts the aligned positions before the next gap: positive
// means a base present only in the reference, negative one present only in the
// query. Aligned stretches come out as 'M' because the delta encoding does not
// say which of them differ; label_operations settles that from the sequence.
std::vector<Op> delta_to_ops(const std::vector<int64_t> &deltas,
                             int64_t ref_span, int64_t qry_span) {
    std::vector<Op> ops;
    int64_t r = 0, q = 0;
    for (const int64_t value : deltas) {
        const int64_t run = std::abs(value) - 1;
        append_op(ops, 'M', run);
        r += run;
        q += run;
        if (value > 0) {
            append_op(ops, 'D', 1);
            r += 1;
        } else {
            append_op(ops, 'I', 1);
            q += 1;
        }
    }
    const int64_t tail_r = ref_span - r, tail_q = qry_span - q;
    if (tail_r != tail_q || tail_r < 0)
        throw std::runtime_error("delta run does not reconcile: " +
                                 std::to_string(tail_r) + " vs " + std::to_string(tail_q));
    append_op(ops, 'M', tail_r);
    return ops;
}

std::map<Key, Alignment> read_delta(const std::string &fname, const std::set<Key> &keys) {
    std::ifstream in(fname);
    if (!in)
        throw std::runtime_error("cannot open " + fname);

    std::map<Key, Alignment> alignments;
    std::string line;
    bool pending = false;
    Key key{};
    char strand = '+';
    std::vector<int64_t> deltas;

    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' ||
                                 line.back() == ' ' || line.back() == '\t'))
            line.pop_back();
        if (line.empty() || line[0] == '>') {
            pending = false;
            continue;
        }

        const std::vector<std::string_view> f = split_ws(line);
        if (f.size() >= 7) {
            // start of an alignment: ref span, query span, error counts
            const int64_t s1 = to_int(f[0]), e1 = to_int(f[1]);
            const int64_t s2 = to_int(f[2]), e2 = to_int(f[3]);
            if (s2 <= e2) {
                strand = '+';
                key = {s1 - 1, e1, s2 - 1, e2};
            } else {
                strand = '-';
                key = {s1 - 1, e1, e2 - 1, s2};
            }
            pending = keys.empty() || keys.find(key) != keys.end();
            deltas.clear();
            continue;
        }
        if (!pending)
            continue;

        const int64_t value = to_int(f[0]);
        if (value == 0) {
            Alignment alignment;
            alignment.strand = strand;
            alignment.has_ops = true;
            alignment.ops = delta_to_ops(deltas, key[1] - key[0], key[3] - key[2]);
            alignments[key] = std::move(alignment);
            pending = false;
        } else {
            deltas.push_back(value);
        }
    }
    return alignments;
}

// Coords give the alignment intervals but no base-level alignment, so blocks
// matched this way still have to take their alignment from SyRI's variants.
std::map<Key, Alignment> read_coords(const std::string &fname, const std::set<Key> &keys) {
    std::ifstream in(fname);
    if (!in)
        throw std::runtime_error("cannot open " + fname);

    std::map<Key, Alignment> alignments;
    std::string line;
    while (std::getline(in, line)) {
        std::replace(line.begin(), line.end(), '|', ' ');
        const std::vector<std::string_view> f = split_ws(line);
        if (f.size() < 4)
            continue;
        int64_t s1, e1, s2, e2;
        try {
            s1 = to_int(f[0]); e1 = to_int(f[1]);
            s2 = to_int(f[2]); e2 = to_int(f[3]);
        } catch (const std::exception &) {
            continue;   // header or separator line
        }
        Alignment alignment;
        Key key;
        if (s2 <= e2) {
            alignment.strand = '+';
            key = {s1 - 1, e1, s2 - 1, e2};
        } else {
            alignment.strand = '-';
            key = {s1 - 1, e1, e2 - 1, s2};
        }
        if (!keys.empty() && keys.find(key) == keys.end())
            continue;
        alignment.has_ops = false;
        alignments[key] = std::move(alignment);
    }
    return alignments;
}

std::string detect_format(const std::string &fname) {
    std::ifstream in(fname);
    if (!in)
        throw std::runtime_error("cannot open " + fname);

    std::string line;
    for (int i = 0; i < 20 && std::getline(in, line); ++i) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();
        if (!line.empty() && line[0] == '>' && split_ws(line).size() == 4)
            return "delta";
        if (line.find("[S1]") != std::string::npos || line.rfind("=====", 0) == 0)
            return "coords";
        const std::vector<std::string_view> f = split(line, '\t');
        if (f.size() >= 12 && (f[4] == "+" || f[4] == "-"))
            return "paf";
    }
    throw std::runtime_error(fname + ": cannot tell whether this is PAF, delta or coords");
}

// ------------------------------------------------------------ inversions ----

std::string invert_query(const std::string &qry_seq, const std::vector<Interval> &regions) {
    std::string out;
    out.reserve(qry_seq.size());
    int64_t prev = 0;
    for (const Interval &r : regions) {
        if (r.begin < prev)
            throw std::runtime_error("inverted regions are not sorted");
        out.append(qry_seq, prev, r.begin - prev);
        for (int64_t i = r.end - 1; i >= r.begin; --i)
            out.push_back(complement(qry_seq[i]));
        prev = r.end;
    }
    out.append(qry_seq, prev, qry_seq.size() - prev);
    return out;
}

std::pair<int64_t, int64_t> remap_query(int64_t qbegin, int64_t qend,
                                        const std::vector<Interval> &regions) {
    for (const Interval &r : regions) {
        if (qbegin >= r.begin && qend <= r.end)
            return {r.begin + r.end - qend, r.begin + r.end - qbegin};
    }
    return {qbegin, qend};
}

// Where a block has to stop so it does not run into a reverse-complemented
// region.
int64_t inverted_limit(int64_t qbegin, int64_t qend, const std::vector<Interval> &regions) {
    int64_t limit = qend;
    for (const Interval &r : regions) {
        if (qbegin < r.begin && r.begin < limit)
            limit = r.begin;
    }
    return limit;
}

// SyRI names alignment records with an 'AL' suffix; those belonging to an
// inverted block start with INV.
bool is_inverted_alignment(const std::string &type) {
    return type.size() >= 2 && type.compare(type.size() - 2, 2, "AL") == 0 &&
           type.rfind("INV", 0) == 0;
}

// --------------------------------------------------------- block walking ----

// Expand one alignment block into operations, applying its variants.
std::vector<Op> block_ops(int64_t rbegin, int64_t rend, int64_t qbegin, int64_t qend,
                          const std::vector<const Variant *> &variants,
                          int64_t *skew, int64_t *dropped) {
    std::vector<Op> ops;
    int64_t r = rbegin, q = qbegin;
    *skew = 0;
    *dropped = 0;

    for (const Variant *v : variants) {
        // INS and DEL are anchored on a shared base that belongs to the
        // surrounding match, so the variant proper starts one base in
        int64_t vr = v->rbegin, vq = v->qbegin;
        if (v->type == VarType::INS || v->type == VarType::DEL) {
            vr += 1;
            vq += 1;
        }
        // variants must advance monotonically; one that starts behind the
        // position already reached overlaps an earlier record
        if (vr < r || vq < q) {
            ++*dropped;
            continue;
        }

        // the stretch before the variant is identical on both sides; if the two
        // sides disagree the block is internally inconsistent, so align what we
        // can and absorb the difference
        const int64_t eq_r = vr - r, eq_q = vq - q;
        append_op(ops, '=', std::min(eq_r, eq_q));
        append_op(ops, 'D', eq_r - eq_q);
        append_op(ops, 'I', eq_q - eq_r);
        *skew += std::abs(eq_r - eq_q);
        r = vr;
        q = vq;

        switch (v->type) {
            case VarType::SNP:
                append_op(ops, 'X', v->rend - v->rbegin);
                r = v->rend;
                q = v->qend;
                break;
            case VarType::DEL:
                append_op(ops, 'D', v->rend - vr);
                r = v->rend;
                break;
            case VarType::INS:
                append_op(ops, 'I', v->qend - vq);
                q = v->qend;
                break;
            case VarType::HDR: {
                // a diverged stretch: pair up what lines up, the rest is indel
                const int64_t rl = v->rend - vr, ql = v->qend - vq;
                append_op(ops, 'X', std::min(rl, ql));
                append_op(ops, 'D', rl - ql);
                append_op(ops, 'I', ql - rl);
                r = v->rend;
                q = v->qend;
                break;
            }
        }
    }

    const int64_t eq_r = rend - r, eq_q = qend - q;
    append_op(ops, '=', std::min(eq_r, eq_q));
    append_op(ops, 'D', eq_r - eq_q);
    append_op(ops, 'I', eq_q - eq_r);
    *skew += std::abs(eq_r - eq_q);
    return ops;
}

// Drop leading operations until the block starts at or after (min_r,min_q), so
// blocks that overlap their predecessor can still be chained.
void trim_front(std::vector<Op> *ops, int64_t rbegin, int64_t qbegin,
                int64_t min_r, int64_t min_q, int64_t *out_r, int64_t *out_q) {
    int64_t r = rbegin, q = qbegin;
    size_t i = 0;
    while (i < ops->size() && (r < min_r || q < min_q)) {
        const char op = (*ops)[i].op;
        const int64_t length = (*ops)[i].length;
        const bool cr = consumes_ref(op), cq = consumes_qry(op);

        int64_t take;
        if ((r < min_r && !cr) || (q < min_q && !cq)) {
            // cannot advance the axis that is still behind, so it all has to go
            take = length;
        } else {
            const int64_t need_r = cr ? std::max<int64_t>(0, min_r - r) : 0;
            const int64_t need_q = cq ? std::max<int64_t>(0, min_q - q) : 0;
            take = std::min(length, std::max(need_r, need_q));
        }

        if (cr) r += take;
        if (cq) q += take;
        if (take >= length) {
            ++i;
        } else {
            (*ops)[i].length = length - take;
            break;
        }
    }
    ops->erase(ops->begin(), ops->begin() + static_cast<std::ptrdiff_t>(i));
    *out_r = r;
    *out_q = q;
}

// Drop trailing operations until the block ends at or before (max_r,max_q).
// Needed where a forward block runs into a region that was reverse-complemented
// for an inversion: those query bases no longer read in the same direction, so
// the overhang cannot stay aligned.
void trim_back(std::vector<Op> *ops, int64_t rend, int64_t qend,
               int64_t max_r, int64_t max_q) {
    int64_t r = rend, q = qend;
    size_t i = ops->size();
    while (i > 0 && (r > max_r || q > max_q)) {
        const char op = (*ops)[i-1].op;
        const int64_t length = (*ops)[i-1].length;
        const bool cr = consumes_ref(op), cq = consumes_qry(op);

        int64_t take;
        if ((r > max_r && !cr) || (q > max_q && !cq)) {
            take = length;
        } else {
            const int64_t need_r = cr ? std::max<int64_t>(0, r - max_r) : 0;
            const int64_t need_q = cq ? std::max<int64_t>(0, q - max_q) : 0;
            take = std::min(length, std::max(need_r, need_q));
        }

        if (cr) r -= take;
        if (cq) q -= take;
        if (take >= length) {
            --i;
        } else {
            (*ops)[i-1].length = length - take;
            break;
        }
    }
    ops->resize(i);
}

// nucmer records where the gaps are but not which aligned bases differ, so that
// distinction has to be recovered from the sequence itself. N positions are
// settled later, by label_operations.
std::vector<Op> resolve_matches(const std::vector<Op> &ops,
                                const std::string &ref_seq, const std::string &qry_seq,
                                int64_t rbegin, int64_t qbegin) {
    bool any = false;
    for (const Op &o : ops)
        any = any || o.op == 'M';
    if (!any)
        return ops;

    std::vector<Op> out;
    int64_t r = rbegin, q = qbegin;
    for (const Op &o : ops) {
        if (o.op != 'M') {
            append_op(out, o.op, o.length);
            if (consumes_ref(o.op)) r += o.length;
            if (consumes_qry(o.op)) q += o.length;
            continue;
        }
        for (int64_t i = 0; i < o.length; ++i)
            append_op(out, ref_seq[r + i] == qry_seq[q + i] ? '=' : 'X', 1);
        r += o.length;
        q += o.length;
    }
    return out;
}

// Align two short stretches whose lengths differ by one unexplained indel.
// Anchors on the longest common prefix and suffix and puts the indel between
// them, which is where a length difference SyRI did not account for actually
// sits. Anything still unmatched in the middle becomes a substitution.
std::vector<Op> realign_window(std::string_view ref, std::string_view qry) {
    std::vector<Op> ops;
    const int64_t n = int64_t(ref.size()), m = int64_t(qry.size());

    int64_t prefix = 0;
    while (prefix < std::min(n, m) && ref[prefix] == qry[prefix])
        ++prefix;

    int64_t suffix = 0;
    while (suffix < std::min(n, m) - prefix && ref[n-1-suffix] == qry[m-1-suffix])
        ++suffix;

    const int64_t mid_r = n - prefix - suffix, mid_q = m - prefix - suffix;
    append_op(ops, '=', prefix);
    append_op(ops, 'X', std::min(mid_r, mid_q));
    append_op(ops, 'D', mid_r - mid_q);
    append_op(ops, 'I', mid_q - mid_r);
    append_op(ops, '=', suffix);
    return ops;
}

// A matching run that does not match means an adjacent indel was placed on the
// wrong side of it, so the run is realigned together with the indels next to
// it. Whatever cannot be recovered that way is emitted as 'X' rather than left
// as a false match.
std::vector<Op> repair_cigar(const std::vector<Op> &cigar,
                             const std::string &ref_seq, const std::string &qry_seq,
                             int64_t *repaired) {
    std::vector<Op> out;
    int64_t r = 0, q = 0;
    size_t i = 0;
    *repaired = 0;

    while (i < cigar.size()) {
        const char op = cigar[i].op;
        const int64_t length = cigar[i].length;

        if (op == '=' && std::memcmp(ref_seq.data() + r, qry_seq.data() + q,
                                     size_t(length)) != 0) {
            int64_t rbegin = r, qbegin = q;
            // an indel just before the run could equally belong just after it
            if (!out.empty() && (out.back().op == 'I' || out.back().op == 'D') &&
                out.back().length <= MAX_SHIFT) {
                if (out.back().op == 'D') rbegin -= out.back().length;
                else                      qbegin -= out.back().length;
                out.pop_back();
            }

            size_t j = i + 1;
            int64_t rend = r + length, qend = q + length;
            while (j < cigar.size() && (cigar[j].op == 'I' || cigar[j].op == 'D') &&
                   cigar[j].length <= MAX_SHIFT) {
                if (cigar[j].op == 'D') rend += cigar[j].length;
                else                    qend += cigar[j].length;
                ++j;
            }

            for (const Op &o : realign_window(
                     std::string_view(ref_seq).substr(rbegin, rend - rbegin),
                     std::string_view(qry_seq).substr(qbegin, qend - qbegin)))
                append_op(out, o.op, o.length);
            ++*repaired;
            r = rend;
            q = qend;
            i = j;
            continue;
        }

        append_op(out, op, length);
        if (consumes_ref(op)) r += length;
        if (consumes_qry(op)) q += length;
        ++i;
    }
    return out;
}

// Relabel every operation that consumes both sequences, so that = and X only
// ever cover positions where neither side has an N, and M covers the rest.
std::vector<Op> label_operations(const std::vector<Op> &cigar,
                                 const std::string &ref_seq, const std::string &qry_seq) {
    std::vector<Op> out;
    int64_t r = 0, q = 0;
    for (const Op &o : cigar) {
        if (o.op != '=' && o.op != 'X' && o.op != 'M') {
            append_op(out, o.op, o.length);
            if (consumes_ref(o.op)) r += o.length;
            if (consumes_qry(o.op)) q += o.length;
            continue;
        }
        for (int64_t i = 0; i < o.length; ++i) {
            const char a = ref_seq[r + i], b = qry_seq[q + i];
            const char label = (a == 'N' || b == 'N') ? 'M' : (a == b ? '=' : 'X');
            append_op(out, label, 1);
        }
        r += o.length;
        q += o.length;
    }
    return out;
}

// --------------------------------------------------------------- stats ------

constexpr int8_t UNCOVERED = 3;

// Per reference position: 0 identical, 1 substituted, 2 absent from the query,
// UNCOVERED where these operations say nothing.
void project_reference(const std::vector<Op> &ops, int64_t obegin, int64_t rbegin,
                       std::vector<int8_t> *state) {
    std::fill(state->begin(), state->end(), UNCOVERED);
    int64_t r = obegin - rbegin;
    for (const Op &o : ops) {
        if (o.op == '=') {
            std::fill(state->begin() + r, state->begin() + r + o.length, int8_t(0));
            r += o.length;
        } else if (o.op == 'X' || o.op == 'M') {
            std::fill(state->begin() + r, state->begin() + r + o.length, int8_t(1));
            r += o.length;
        } else if (o.op == 'D') {
            std::fill(state->begin() + r, state->begin() + r + o.length, int8_t(2));
            r += o.length;
        }
    }
}

struct Stats {
    int64_t blocks = 0, used = 0, skipped = 0;
    int64_t skew = 0, dropped = 0;
    int64_t inverted = 0, inverted_bp = 0, merges = 0;
    int64_t clipped_by_inversion = 0;
    int64_t repaired = 0;
    int64_t from_aln = 0, from_syri = 0, matched_intervals = 0;
    int64_t eq_divergence = 0;
    int64_t syri_unsupported = 0, syri_contradicted = 0;
    std::string aln_format;
    std::string qry_header;
};

}  // namespace

// ----------------------------------------------------------------- build ----

namespace {

struct Result {
    std::vector<Op> cigar;
    std::string ref_seq;
    std::string qry_seq;
    Stats stats;
};

Result build_cigar(const std::string &ref_fname, const std::string &qry_fname,
                   const std::string &syri_fname, const std::string &aln_fname,
                   std::string aln_format) {
    Result result;
    std::string ref_header;
    result.ref_seq = read_fasta(ref_fname, &ref_header);
    std::string qry_seq = read_fasta(qry_fname, &result.stats.qry_header);

    std::unordered_map<std::string, int> parent_ids;
    auto intern = [&parent_ids](std::string_view name) {
        const std::string key(name);
        auto it = parent_ids.find(key);
        if (it != parent_ids.end())
            return it->second;
        const int id = int(parent_ids.size());
        parent_ids.emplace(key, id);
        return id;
    };

    std::vector<Block> blocks;
    std::vector<Variant> variants;
    std::vector<Interval> inverted;

    stream_syri(syri_fname, [&](int64_t rbegin, int64_t rend, int64_t qbegin, int64_t qend,
                                std::string_view id, std::string_view parent,
                                std::string_view type) {
        // inverted records are written with the query interval running backwards
        if (qend > 0 && qbegin + 1 > qend) {
            const int64_t b = qend - 1, e = qbegin + 1;
            qbegin = b;
            qend = e;
        }
        const std::string type_str(type);
        // flip exactly the intervals whose alignments we intend to chain, so the
        // region reverse-complemented is the one the CIGAR describes
        if (INVERTED_TYPES.count(type_str))
            inverted.push_back({qbegin, qend});
        if (CHAIN_TYPES.count(type_str)) {
            blocks.push_back({rbegin, rend, qbegin, qend, std::string(id), type_str,
                              intern(parent), Key{rbegin, rend, qbegin, qend}});
        } else if (type == "SNP" || type == "INS" || type == "DEL" || type == "HDR") {
            const VarType vt = type == "SNP" ? VarType::SNP
                             : type == "INS" ? VarType::INS
                             : type == "DEL" ? VarType::DEL
                                             : VarType::HDR;
            variants.push_back({rbegin, rend, qbegin, qend, vt, intern(parent)});
        }
    });

    // Overlapping inversions cannot each be flipped, so merge them. The merged
    // span still reverse-complements every alignment inside it correctly; what
    // is lost is the order of the two inversions relative to each other, which
    // trimming resolves during chaining.
    std::sort(inverted.begin(), inverted.end(),
              [](const Interval &a, const Interval &b) { return a.begin < b.begin; });
    std::vector<Interval> merged;
    for (const Interval &r : inverted) {
        if (!merged.empty() && r.begin < merged.back().end) {
            merged.back().end = std::max(merged.back().end, r.end);
            ++result.stats.merges;
        } else {
            merged.push_back(r);
        }
    }
    inverted = merged;
    if (result.stats.merges)
        std::cerr << "WARNING: merged " << result.stats.merges
                  << " overlapping inverted regions" << std::endl;

    result.qry_seq = invert_query(qry_seq, inverted);
    qry_seq.clear();
    qry_seq.shrink_to_fit();
    const std::string &ref_seq = result.ref_seq;
    const std::string &qry = result.qry_seq;
    const int64_t ref_len = int64_t(ref_seq.size()), qry_len = int64_t(qry.size());

    result.stats.inverted = int64_t(inverted.size());
    for (const Interval &r : inverted)
        result.stats.inverted_bp += r.end - r.begin;

    for (Block &b : blocks) {
        const auto [nb, ne] = remap_query(b.qbegin, b.qend, inverted);
        b.qbegin = nb;
        b.qend = ne;
    }
    for (Variant &v : variants) {
        const auto [nb, ne] = remap_query(v.qbegin, v.qend, inverted);
        v.qbegin = nb;
        v.qend = ne;
    }

    std::sort(variants.begin(), variants.end(), [](const Variant &a, const Variant &b) {
        return std::tie(a.rbegin, a.rend, a.qbegin, a.qend) <
               std::tie(b.rbegin, b.rend, b.qbegin, b.qend);
    });
    std::sort(blocks.begin(), blocks.end(), [](const Block &a, const Block &b) {
        return std::tie(a.rbegin, a.rend, a.qbegin, a.qend) <
               std::tie(b.rbegin, b.rend, b.qbegin, b.qend);
    });
    result.stats.blocks = int64_t(blocks.size());

    // Read the alignments only now, so only the blocks that will actually be
    // chained get their CIGARs parsed. A whole-genome PAF can be hundreds of
    // megabytes of CIGAR for thousands of records when a few hundred are wanted.
    std::map<Key, Alignment> alignments;
    if (!aln_fname.empty()) {
        if (aln_format.empty())
            aln_format = detect_format(aln_fname);
        std::set<Key> keys;
        for (const Block &b : blocks)
            keys.insert(b.key);
        if (aln_format == "paf")        alignments = read_paf(aln_fname, keys);
        else if (aln_format == "delta") alignments = read_delta(aln_fname, keys);
        else                            alignments = read_coords(aln_fname, keys);
    }
    result.stats.aln_format = aln_format;

    // A variant and the alignment it belongs to are both children of the same
    // SyRI block, so match them on the parent id first. Going by coordinates
    // alone would hand a variant to whichever block happens to span it, which
    // neighbouring blocks routinely do where they overlap.
    std::unordered_map<int, std::vector<const Block *>> by_parent;
    for (const Block &b : blocks)
        by_parent[b.parent].push_back(&b);

    std::unordered_map<std::string, std::vector<const Variant *>> per_block;
    for (const Variant &v : variants) {
        auto it = by_parent.find(v.parent);
        if (it == by_parent.end())
            continue;
        for (const Block *b : it->second) {
            if (v.rbegin >= b->rbegin && v.rend <= b->rend &&
                v.qbegin >= b->qbegin && v.qend <= b->qend) {
                per_block[b->id].push_back(&v);
                break;
            }
        }
    }

    std::vector<Op> cigar;
    int64_t r = 0, q = 0;
    const std::vector<const Variant *> no_variants;
    std::vector<int8_t> state_a, state_s;

    for (const Block &b : blocks) {
        auto it = per_block.find(b.id);
        const std::vector<const Variant *> &block_variants =
            it == per_block.end() ? no_variants : it->second;

        int64_t skew = 0, dropped = 0;
        std::vector<Op> syri_ops = block_ops(b.rbegin, b.rend, b.qbegin, b.qend,
                                             block_variants, &skew, &dropped);

        auto entry = alignments.find(b.key);
        const bool have_alignment = entry != alignments.end() && entry->second.has_ops;

        std::vector<Op> ops;
        if (have_alignment) {
            // the aligner's own output carries the real alignment, so use it
            // rather than reconstructing one from SyRI's variant records
            ops = entry->second.ops;
            skew = dropped = 0;
            ++result.stats.from_aln;

            const auto [span_r, span_q] = ops_span(ops);
            if (span_r != b.rend - b.rbegin || span_q != b.qend - b.qbegin)
                std::cerr << "WARNING: alignment for " << b.id << " spans " << span_r
                          << "/" << span_q << ", block is " << (b.rend - b.rbegin)
                          << "/" << (b.qend - b.qbegin) << std::endl;
        } else {
            if (entry != alignments.end())
                ++result.stats.matched_intervals;
            ops = syri_ops;
            ++result.stats.from_syri;
        }

        // trim both the alignment and the SyRI reconstruction to the same window
        // so the cross-check below compares like with like
        int64_t brbegin, bqbegin, srbegin, sqbegin;
        trim_front(&ops, b.rbegin, b.qbegin, r, q, &brbegin, &bqbegin);
        trim_front(&syri_ops, b.rbegin, b.qbegin, r, q, &srbegin, &sqbegin);
        (void)sqbegin;

        if (!ops.empty() && !is_inverted_alignment(b.type)) {
            // A forward block must not reach into a region flipped for an
            // inversion; the inverted alignment itself owns that sequence. Note
            // that a palindromic repeat can align in both orientations, so the
            // sequence check cannot tell the two readings apart here -- the
            // clipped total below is how much sequence that choice affects.
            const int64_t limit = inverted_limit(bqbegin, b.qend, inverted);
            if (limit < b.qend) {
                const int64_t before = ops_span(ops).second;
                trim_back(&ops, b.rend, b.qend, ref_len, limit);
                trim_back(&syri_ops, b.rend, b.qend, ref_len, limit);
                result.stats.clipped_by_inversion += before - ops_span(ops).second;
            }
        }
        if (ops.empty()) {
            ++result.stats.skipped;
            continue;
        }
        result.stats.skew += skew;
        result.stats.dropped += dropped;
        ++result.stats.used;

        // only now that the surviving window is known is it worth resolving
        // which aligned bases match; doing it earlier would compare against
        // sequence that is about to be discarded, possibly reverse-complemented
        ops = resolve_matches(ops, ref_seq, qry, brbegin, bqbegin);

        if (have_alignment) {
            int64_t aln_eq = 0, syri_eq = 0;
            for (const Op &o : ops)      if (o.op == '=') aln_eq += o.length;
            for (const Op &o : syri_ops) if (o.op == '=') syri_eq += o.length;
            result.stats.eq_divergence += std::abs(aln_eq - syri_eq);

            const int64_t span = b.rend - b.rbegin;
            state_a.assign(size_t(span), UNCOVERED);
            state_s.assign(size_t(span), UNCOVERED);
            project_reference(ops, brbegin, b.rbegin, &state_a);
            project_reference(syri_ops, srbegin, b.rbegin, &state_s);
            for (int64_t i = 0; i < span; ++i) {
                if (state_a[i] == UNCOVERED || state_s[i] == UNCOVERED)
                    continue;
                if (state_s[i] == 0 && state_a[i] != 0) ++result.stats.syri_unsupported;
                if (state_a[i] == 0 && state_s[i] != 0) ++result.stats.syri_contradicted;
            }
        }

        // unalignable sequence between the previous block and this one
        append_op(cigar, 'D', brbegin - r);
        append_op(cigar, 'I', bqbegin - q);
        r = brbegin;
        q = bqbegin;
        for (const Op &o : ops) {
            append_op(cigar, o.op, o.length);
            if (consumes_ref(o.op)) r += o.length;
            if (consumes_qry(o.op)) q += o.length;
        }
    }

    append_op(cigar, 'D', ref_len - r);
    append_op(cigar, 'I', qry_len - q);

    cigar = repair_cigar(cigar, ref_seq, qry, &result.stats.repaired);
    // finally apply the N convention, whatever the operations came from
    result.cigar = label_operations(cigar, ref_seq, qry);
    return result;
}

// Check the CIGAR is global and that every label is justified by the sequences.
struct Verification {
    int64_t ref_consumed = 0, qry_consumed = 0;
    int64_t bad_eq = 0, bad_x = 0, bad_m = 0;
};

Verification verify(const std::vector<Op> &cigar,
                    const std::string &ref_seq, const std::string &qry_seq) {
    Verification v;
    int64_t r = 0, q = 0;
    for (const Op &o : cigar) {
        if (o.op == '=' || o.op == 'X' || o.op == 'M') {
            for (int64_t i = 0; i < o.length; ++i) {
                const char a = ref_seq[r + i], b = qry_seq[q + i];
                const bool has_N = (a == 'N' || b == 'N');
                if (o.op == '=' && (has_N || a != b)) ++v.bad_eq;
                if (o.op == 'X' && (has_N || a == b)) ++v.bad_x;
                if (o.op == 'M' && !has_N)            ++v.bad_m;
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
    "usage: syri_to_cigar [--format FMT] REF.fa QRY.fa SYRI.out ALIGNMENTS OUT_PREFIX\n"
    "\n"
    "Build a global CIGAR alignment between a reference and a query from SyRI output,\n"
    "taking the base-level alignment from the alignments SyRI was run on.\n"
    "\n"
    "positional arguments:\n"
    "  REF.fa       reference FASTA, single sequence\n"
    "  QRY.fa       query FASTA, single sequence\n"
    "  SYRI.out     SyRI annotation table, used to classify each alignment\n"
    "  ALIGNMENTS   the alignments SyRI was run on, in one of the formats below\n"
    "  OUT_PREFIX   prefix for the output files, including any directory\n"
    "\n"
    "options:\n"
    "  --format FMT  paf, delta or coords. Detected from the file if not given.\n"
    "  -h, --help    show this message\n"
    "\n"
    "alignment formats:\n"
    "  paf     minimap2 PAF. Needs CIGARs in the cg:Z: tag, which means --eqx or\n"
    "          -c; records without one fall back to the SyRI variant records.\n"
    "  delta   nucmer .delta. Carries gap positions but not which aligned bases\n"
    "          differ, so = and X are resolved against the sequences.\n"
    "  coords  show-coords output, padded or -T. Gives only intervals, so every\n"
    "          block still takes its alignment from the SyRI variant records; use\n"
    "          delta instead if you have it.\n"
    "\n"
    "output:\n"
    "  OUT_PREFIX.cigar   the CIGAR, as run-length pairs such as 135=5I6=1X\n"
    "  OUT_PREFIX.qry.fa  the query with inverted regions reverse-complemented. The\n"
    "                     CIGAR aligns the reference against THIS sequence, not the\n"
    "                     original query, so validating against QRY.fa will fail.\n"
    "  A summary of the chaining, and checks that the alignment is global and that\n"
    "  every =, X and M is correct, are written to standard error.\n"
    "\n"
    "operations:\n"
    "  =  identical bases, neither an N, consumes reference and query\n"
    "  X  differing bases, neither an N, consumes reference and query\n"
    "  M  a position where either side is an N, consumes reference and query\n"
    "  D  consumes reference only\n"
    "  I  consumes query only\n"
    "\n"
    "examples:\n"
    "  ./syri_to_cigar ref.fa qry.fa syri.out aln.paf out/chr22\n"
    "  ./syri_to_cigar ref.fa qry.fa syri.out out.delta out/chr22\n"
    "  writes out/chr22.cigar and out/chr22.qry.fa\n";

}  // namespace

int main(int argc, char **argv) {
    std::vector<std::string> args;
    std::string aln_format;
    bool want_help = false;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help") {
            want_help = true;
        } else if (a == "--format" && i + 1 < argc) {
            aln_format = argv[++i];
        } else if (a.rfind("--format=", 0) == 0) {
            aln_format = a.substr(9);
        } else {
            args.push_back(a);
        }
    }

    if (want_help || args.size() != 5) {
        std::fputs(USAGE, stderr);
        return want_help ? 0 : 2;
    }
    if (!aln_format.empty() && aln_format != "paf" && aln_format != "delta" &&
        aln_format != "coords") {
        std::cerr << "unknown format '" << aln_format
                  << "'; expected one of paf, delta, coords" << std::endl;
        return 2;
    }

    try {
        const Result result = build_cigar(args[0], args[1], args[2], args[3], aln_format);
        const Stats &s = result.stats;

        std::map<char, int64_t> counts;
        for (const Op &o : result.cigar)
            counts[o.op] += o.length;

        const Verification v = verify(result.cigar, result.ref_seq, result.qry_seq);
        const int64_t ref_len = int64_t(result.ref_seq.size());
        const int64_t qry_len = int64_t(result.qry_seq.size());

        std::cerr << "inverted regions flipped: " << s.inverted
                  << " (" << s.inverted_bp << " bp)\n";
        if (s.merges)
            std::cerr << "  overlapping inverted regions merged: " << s.merges << "\n";
        if (s.clipped_by_inversion)
            std::cerr << "  query bases clipped from forward blocks reaching into them: "
                      << s.clipped_by_inversion << "\n";
        std::cerr << "blocks chained: " << s.used << " / " << s.blocks
                  << " (dropped as fully overlapping: " << s.skipped << ")\n"
                  << "  alignment from " << (s.aln_format.empty() ? "none" : s.aln_format)
                  << ": " << s.from_aln << ", reconstructed from SyRI variants: "
                  << s.from_syri << "\n"
                  << "  '=' disagreement between alignment and SyRI: " << s.eq_divergence << "\n"
                  << "    SyRI-implied matches with no CIGAR support: "
                  << s.syri_unsupported << "\n"
                  << "    CIGAR matches contradicted by a SyRI record: "
                  << s.syri_contradicted << "\n";
        if (s.matched_intervals)
            std::cerr << "  intervals confirmed but carrying no alignment: "
                      << s.matched_intervals << "\n";
        std::cerr << "variants dropped as overlapping: " << s.dropped << "\n"
                  << "internal skew reconciled: " << s.skew << "\n"
                  << "runs realigned by repair pass: " << s.repaired << "\n"
                  << "operations: " << result.cigar.size() << "\n";
        for (const char op : std::string("=XMID"))
            std::cerr << "  " << op << ": " << counts[op] << "\n";
        std::cerr << "reference consumed: " << v.ref_consumed << " / " << ref_len
                  << (v.ref_consumed == ref_len ? " OK" : " MISMATCH") << "\n"
                  << "query consumed:     " << v.qry_consumed << " / " << qry_len
                  << (v.qry_consumed == qry_len ? " OK" : " MISMATCH") << "\n"
                  << "'=' positions that are not identical or involve an N: " << v.bad_eq << "\n"
                  << "'X' positions that are identical or involve an N:     " << v.bad_x << "\n"
                  << "'M' positions with no N on either side:               " << v.bad_m << "\n";
        const int64_t minlen = std::min(ref_len, qry_len);
        std::fprintf(stderr, "identity: %lld / %lld (%.2f%%)\n",
                     (long long)counts['='], (long long)minlen,
                     100.0 * double(counts['=']) / double(minlen));

        const std::string cigar_fname = args[4] + ".cigar";
        const std::string flipped_fname = args[4] + ".qry.fa";

        {
            std::ofstream out(cigar_fname);
            if (!out)
                throw std::runtime_error("cannot write " + cigar_fname);
            for (const Op &o : result.cigar)
                out << o.length << o.op;
            out << "\n";
        }
        {
            // the CIGAR aligns against the query with its inversions flipped, so
            // that sequence has to travel with it
            std::ofstream out(flipped_fname);
            if (!out)
                throw std::runtime_error("cannot write " + flipped_fname);
            out << ">inverted_" << s.qry_header << "\n";
            for (size_t i = 0; i < result.qry_seq.size(); i += 60)
                out << result.qry_seq.substr(i, 60) << "\n";
        }
        std::cerr << "wrote " << cigar_fname << " and " << flipped_fname << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
