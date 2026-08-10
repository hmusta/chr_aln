#!/usr/bin/env python3
"""Build a global CIGAR alignment between a reference and a query from SyRI output.

Usage:
    ./syri_to_cigar.py REF.fa QRY.fa SYRI.out ALIGNMENTS.paf [OUT.cigar]

SyRI classifies each alignment (syntenic, inverted, duplicated, translocated),
which is what decides whether a block can join a colinear chain, but it only
summarises the alignment itself as variant records. The PAF it was run on
carries the real base-level CIGAR, so the classification is taken from SyRI and
the alignment from the PAF. Blocks with no matching PAF record fall back to
reconstructing the alignment from the variant records.

Operations, in the usual convention:
    =   identical bases, consumes reference and query
    X   substituted bases, consumes reference and query
    D   consumes reference only
    I   consumes query only

The alignment is global: the operations consuming reference sum to the reference
length and those consuming query sum to the query length.

Inverted regions are reverse-complemented in the query before chaining, which
makes them colinear with the reference so they contribute matches rather than
indels. Every query coordinate inside such a region is remapped to match. The
reverse-complemented query is written alongside the CIGAR when an output file is
given, since the CIGAR describes the alignment against that sequence.

A CIGAR is still colinear, so translocations and duplicate copies cannot be
reached and come out as D/I, as does anything SyRI left unaligned.
"""

import re
import sys

from helpers import _COMPLEMENT,read_fasta
from check_syri import stream_syri

# blocks that can carry matches once inversions have been flipped
CHAIN_TYPES = ("SYNAL","INVAL")
# regions reverse-complemented in the query. Inverted duplications are left
# alone: they are surplus copies whose reference is already spoken for, so
# flipping them would not make them chainable.
INVERTED_TYPES = ("INVAL","INVTRAL")
VARIANT_TYPES = ("SNP","INS","DEL","HDR")

def read_paf(fname):
    """Index the PAF alignments by their coordinates, keeping the CIGAR.

    The PAF is the alignment SyRI was run on, so its records carry the exact
    base-level alignment that SyRI only summarises as variant records. Query
    coordinates are always on the forward strand; for a '-' record the CIGAR
    describes the reverse complement of that interval, which is precisely the
    sequence we get after flipping an inverted region.
    """
    alignments = {}
    for line in open(fname):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12:
            continue
        qbegin,qend = int(fields[2]),int(fields[3])
        strand = fields[4]
        rbegin,rend = int(fields[7]),int(fields[8])

        cigar = None
        for tag in fields[12:]:
            if tag.startswith("cg:Z:"):
                cigar = tag[5:]
                break
        if cigar is None:
            continue

        ops = [[op,int(n)] for n,op in re.findall(r"(\d+)([=XID])",cigar)]
        alignments[(rbegin,rend,qbegin,qend)] = (strand,ops)
    return alignments

def ops_span(ops):
    return (sum(l for o,l in ops if o in "=XD"),
            sum(l for o,l in ops if o in "=XI"))

def append_op(ops, op, length):
    if length <= 0:
        return
    if ops and ops[-1][0] == op:
        ops[-1][1] += length
    else:
        ops.append([op,length])

def invert_query(qry_seq, regions):
    """Reverse-complement each region of the query in place (by rebuilding it)."""
    pieces = []
    prev = 0
    for begin,end in regions:
        assert(begin >= prev)
        pieces.append(qry_seq[prev:begin])
        pieces.append(qry_seq[begin:end].translate(_COMPLEMENT)[::-1])
        prev = end
    pieces.append(qry_seq[prev:])
    return "".join(pieces)

def remap_query(qbegin, qend, regions):
    """Map a query interval through the reverse-complementing of its region."""
    for begin,end in regions:
        if qbegin >= begin and qend <= end:
            return begin + end - qend,begin + end - qbegin
    return qbegin,qend

def block_ops(rbegin, rend, qbegin, qend, variants):
    """Expand one alignment block into operations, applying its variants.

    Returns the operation list, the number of positions that had to be
    reconciled because the reference and query ran out of step, and the number
    of variants dropped for overlapping one already emitted.
    """
    ops = []
    r,q = rbegin,qbegin
    skew = 0
    dropped = 0

    for vrbegin,vrend,vqbegin,vqend,vtype,*_ in variants:
        # INS and DEL are anchored on a shared base that belongs to the
        # surrounding match, so the variant proper starts one base in
        if vtype in ("INS","DEL"):
            vr,vq = vrbegin + 1,vqbegin + 1
        else:
            vr,vq = vrbegin,vqbegin

        # variants must advance monotonically; one that starts behind the
        # position already reached overlaps an earlier record
        if vr < r or vq < q:
            dropped += 1
            continue

        # the stretch before the variant is identical on both sides; if the two
        # sides disagree the block is internally inconsistent, so align what we
        # can and absorb the difference
        eq_r,eq_q = vr - r,vq - q
        append_op(ops,'=',min(eq_r,eq_q))
        append_op(ops,'D',eq_r - eq_q)
        append_op(ops,'I',eq_q - eq_r)
        skew += abs(eq_r - eq_q)
        r,q = vr,vq

        if vtype == "SNP":
            append_op(ops,'X',vrend - vrbegin)
            r,q = vrend,vqend
        elif vtype == "DEL":
            append_op(ops,'D',vrend - vr)
            r = vrend
        elif vtype == "INS":
            append_op(ops,'I',vqend - vq)
            q = vqend
        elif vtype == "HDR":
            # a diverged stretch: pair up what lines up, the rest is indel
            rl,ql = vrend - vr,vqend - vq
            append_op(ops,'X',min(rl,ql))
            append_op(ops,'D',rl - ql)
            append_op(ops,'I',ql - rl)
            r,q = vrend,vqend

    eq_r,eq_q = rend - r,qend - q
    append_op(ops,'=',min(eq_r,eq_q))
    append_op(ops,'D',eq_r - eq_q)
    append_op(ops,'I',eq_q - eq_r)
    skew += abs(eq_r - eq_q)

    return ops,skew,dropped

def trim_front(ops, rbegin, qbegin, min_r, min_q):
    """Drop leading operations until the block starts at or after (min_r,min_q),
    so blocks that overlap their predecessor can still be chained."""
    r,q = rbegin,qbegin
    i = 0
    while i < len(ops) and (r < min_r or q < min_q):
        op,length = ops[i]
        consumes_r = op in "=XD"
        consumes_q = op in "=XI"

        if (r < min_r and not consumes_r) or (q < min_q and not consumes_q):
            # cannot advance the axis that is still behind, so it all has to go
            take = length
        else:
            need_r = max(0,min_r - r) if consumes_r else 0
            need_q = max(0,min_q - q) if consumes_q else 0
            take = min(length,max(need_r,need_q))

        r += take if consumes_r else 0
        q += take if consumes_q else 0
        if take >= length:
            i += 1
        else:
            ops[i] = [op,length - take]
            break

    return ops[i:],r,q

def trim_back(ops, rend, qend, max_r, max_q):
    """Drop trailing operations until the block ends at or before (max_r,max_q).

    Needed where a forward block runs into a region that was reverse-complemented
    for an inversion: those query bases no longer read in the same direction, so
    the overhang cannot stay aligned.
    """
    r,q = rend,qend
    i = len(ops)
    while i > 0 and (r > max_r or q > max_q):
        op,length = ops[i-1]
        consumes_r = op in "=XD"
        consumes_q = op in "=XI"

        if (r > max_r and not consumes_r) or (q > max_q and not consumes_q):
            take = length
        else:
            need_r = max(0,r - max_r) if consumes_r else 0
            need_q = max(0,q - max_q) if consumes_q else 0
            take = min(length,max(need_r,need_q))

        r -= take if consumes_r else 0
        q -= take if consumes_q else 0
        if take >= length:
            i -= 1
        else:
            ops[i-1] = [op,length - take]
            break

    return ops[:i],r,q

def inverted_limit(qbegin, qend, regions):
    """Where a block has to stop so it does not run into a reverse-complemented
    region."""
    limit = qend
    for begin,end in regions:
        if qbegin < begin < limit:
            limit = begin
    return limit

def realign_window(ref, qry):
    """Align two short stretches whose lengths differ by one unexplained indel.

    Anchors on the longest common prefix and suffix and puts the indel between
    them, which is where a length difference SyRI did not account for actually
    sits. Anything still unmatched in the middle becomes a substitution.
    """
    ops = []
    n,m = len(ref),len(qry)

    prefix = 0
    while prefix < min(n,m) and ref[prefix] == qry[prefix]:
        prefix += 1

    suffix = 0
    while suffix < min(n,m) - prefix and ref[n-1-suffix] == qry[m-1-suffix]:
        suffix += 1

    mid_r,mid_q = n - prefix - suffix,m - prefix - suffix
    append_op(ops,'=',prefix)
    append_op(ops,'X',min(mid_r,mid_q))
    append_op(ops,'D',mid_r - mid_q)
    append_op(ops,'I',mid_q - mid_r)
    append_op(ops,'=',suffix)
    return ops

def repair_cigar(cigar, ref_seq, qry_seq, max_shift=4096):
    """Make every '=' operation honest.

    A matching run that does not match means an adjacent indel was placed on the
    wrong side of it, so the run is realigned together with the indels next to
    it. Whatever cannot be recovered that way is emitted as 'X' rather than left
    as a false match.
    """
    out = []
    r = q = 0
    i = 0
    repaired = 0

    while i < len(cigar):
        op,length = cigar[i]

        if op == '=' and ref_seq[r:r+length] != qry_seq[q:q+length]:
            rbegin,qbegin = r,q
            # an indel just before the run could equally belong just after it
            if out and out[-1][0] in "ID" and out[-1][1] <= max_shift:
                prev_op,prev_len = out.pop()
                rbegin -= prev_len if prev_op == 'D' else 0
                qbegin -= prev_len if prev_op == 'I' else 0

            j = i + 1
            rend,qend = r + length,q + length
            while j < len(cigar) and cigar[j][0] in "ID" and cigar[j][1] <= max_shift:
                rend += cigar[j][1] if cigar[j][0] == 'D' else 0
                qend += cigar[j][1] if cigar[j][0] == 'I' else 0
                j += 1

            for o,l in realign_window(ref_seq[rbegin:rend],qry_seq[qbegin:qend]):
                append_op(out,o,l)
            repaired += 1
            r,q = rend,qend
            i = j
            continue

        append_op(out,op,length)
        r += length if op in "=XD" else 0
        q += length if op in "=XI" else 0
        i += 1

    return out,repaired

def build_cigar(ref_fname, qry_fname, syri_fname, paf_fname=None):
    ref_header,ref_seq = read_fasta(ref_fname)
    qry_header,qry_seq = read_fasta(qry_fname)
    paf = read_paf(paf_fname) if paf_fname else {}

    records = []
    inverted = []
    for entry in stream_syri(syri_fname):
        i,rh,rbegin,rend,rseq,qh,qbegin,qend,qseq,rid,pid,atype,in_inv,copy_type = entry
        if qend > 0 and qbegin + 1 > qend:
            qbegin,qend = qend - 1,qbegin + 1
        records.append((rbegin,rend,qbegin,qend,rid,pid,atype))
        # flip exactly the intervals whose alignments we intend to chain, so the
        # region reverse-complemented is the one the CIGAR describes
        if atype in INVERTED_TYPES:
            inverted.append((qbegin,qend))

    inverted.sort()
    for a,b in zip(inverted,inverted[1:]):
        assert(a[1] <= b[0]),"inverted regions overlap"

    qry_seq = invert_query(qry_seq,inverted)
    ref_len,qry_len = len(ref_seq),len(qry_seq)

    blocks = []
    variants = []
    for rbegin,rend,qbegin,qend,rid,pid,atype in records:
        # the PAF is keyed on the original, unflipped query coordinates
        paf_key = (rbegin,rend,qbegin,qend)
        qbegin,qend = remap_query(qbegin,qend,inverted)
        if atype in CHAIN_TYPES:
            blocks.append((rbegin,rend,qbegin,qend,rid,pid,atype,paf_key))
        elif atype in VARIANT_TYPES:
            variants.append((rbegin,rend,qbegin,qend,atype,pid))

    variants.sort()
    blocks.sort()

    # A variant and the alignment it belongs to are both children of the same
    # SyRI block, so match them on the parent id first. Going by coordinates
    # alone would hand a variant to whichever block happens to span it, which
    # neighbouring blocks routinely do where they overlap.
    by_parent = {}
    for b in blocks:
        by_parent.setdefault(b[5],[]).append(b)

    per_block = {}
    for v in variants:
        for b in by_parent.get(v[5],()):
            if v[0] >= b[0] and v[1] <= b[1] and v[2] >= b[2] and v[3] <= b[3]:
                per_block.setdefault(b[4],[]).append(v)
                break

    cigar = []
    r,q = 0,0
    used = skipped = 0
    total_skew = total_dropped = 0
    from_paf = from_syri = 0
    eq_divergence = 0

    for rbegin,rend,qbegin,qend,rid,pid,atype,paf_key in blocks:
        entry = paf.get(paf_key)
        if entry is not None:
            # the PAF carries the real alignment, so use it rather than
            # reconstructing one from SyRI's variant records
            strand,paf_ops = entry
            ops = [[op,length] for op,length in paf_ops]
            skew = dropped = 0
            from_paf += 1

            span_r,span_q = ops_span(ops)
            if (span_r,span_q) != (rend - rbegin,qend - qbegin):
                print(f"WARNING: PAF CIGAR for {rid} spans {span_r}/{span_q}, "
                      f"block is {rend-rbegin}/{qend-qbegin}",file=sys.stderr)
            # cross-check what SyRI's variants alone would have produced
            syri_ops,_,_ = block_ops(rbegin,rend,qbegin,qend,per_block.get(rid,[]))
            paf_eq = sum(l for o,l in ops if o == '=')
            syri_eq = sum(l for o,l in syri_ops if o == '=')
            eq_divergence += abs(paf_eq - syri_eq)
        else:
            ops,skew,dropped = block_ops(rbegin,rend,qbegin,qend,per_block.get(rid,[]))
            from_syri += 1

        ops,brbegin,bqbegin = trim_front(ops,rbegin,qbegin,r,q)
        if ops and atype not in ("INVAL","INVTRAL"):
            # a forward block must not reach into a region flipped for an
            # inversion; the inverted alignment itself owns that sequence
            limit = inverted_limit(bqbegin,qend,inverted)
            if limit < qend:
                ops,_,_ = trim_back(ops,rend,qend,ref_len,limit)
        if not ops:
            skipped += 1
            continue
        total_skew += skew
        total_dropped += dropped
        used += 1

        # unalignable sequence between the previous block and this one
        append_op(cigar,'D',brbegin - r)
        append_op(cigar,'I',bqbegin - q)
        r,q = brbegin,bqbegin
        for op,length in ops:
            append_op(cigar,op,length)
            r += length if op in "=XD" else 0
            q += length if op in "=XI" else 0

    append_op(cigar,'D',ref_len - r)
    append_op(cigar,'I',qry_len - q)

    cigar,repaired = repair_cigar(cigar,ref_seq,qry_seq)

    stats = dict(blocks=len(blocks),used=used,skipped=skipped,
                 skew=total_skew,dropped=total_dropped,inverted=len(inverted),
                 inverted_bp=sum(e - b for b,e in inverted),repaired=repaired,
                 from_paf=from_paf,from_syri=from_syri,eq_divergence=eq_divergence)
    return cigar,ref_seq,qry_seq,stats

def verify(cigar, ref_seq, qry_seq):
    """Check the CIGAR is global and that every '=' really is identical."""
    r = q = 0
    bad_eq = bad_x = 0
    for op,length in cigar:
        if op == '=':
            if ref_seq[r:r+length] != qry_seq[q:q+length]:
                bad_eq += sum(1 for a,b in zip(ref_seq[r:r+length],qry_seq[q:q+length]) if a != b)
            r += length; q += length
        elif op == 'X':
            bad_x += sum(1 for a,b in zip(ref_seq[r:r+length],qry_seq[q:q+length]) if a == b)
            r += length; q += length
        elif op == 'D':
            r += length
        elif op == 'I':
            q += length
    return r,q,bad_eq,bad_x

USAGE = """usage: syri_to_cigar.py REF.fa QRY.fa SYRI.out ALIGNMENTS.paf [OUT.cigar]

Build a global CIGAR alignment between a reference and a query from SyRI output,
taking the base-level alignment from the PAF that SyRI was run on.

positional arguments:
  REF.fa           reference FASTA, single sequence
  QRY.fa           query FASTA, single sequence
  SYRI.out         SyRI annotation table, used to classify each alignment
  ALIGNMENTS.paf   the PAF SyRI was run on; must carry CIGARs in the cg:Z: tag
  OUT.cigar        where to write the CIGAR (default: standard output)

output:
  OUT.cigar        the CIGAR, as run-length pairs such as 135=5I6=1X
  OUT.cigar.qry.fa the query with inverted regions reverse-complemented. The
                   CIGAR aligns the reference against THIS sequence, not the
                   original query, so validating against QRY.fa will fail.
  A summary of the chaining, and checks that the alignment is global and that
  every = and X is correct, are written to standard error.

operations:
  =  identical bases, consumes reference and query
  X  substituted bases, consumes reference and query
  D  consumes reference only
  I  consumes query only

example:
  ./syri_to_cigar.py ref.fa qry.fa syri.out aln.paf out.cigar
"""

if __name__ == "__main__":
    args = sys.argv[1:]
    if any(a in ("-h","--help") for a in args) or not 4 <= len(args) <= 5:
        sys.stderr.write(USAGE)
        sys.exit(0 if args and args[0] in ("-h","--help") else 2)

    ref_fname,qry_fname,syri_fname,paf_fname = args[:4]
    out_fname = args[4] if len(args) > 4 else None

    cigar,ref_seq,qry_seq,stats = build_cigar(ref_fname,qry_fname,syri_fname,paf_fname)

    counts = {}
    for op,length in cigar:
        counts[op] = counts.get(op,0) + length

    r,q,bad_eq,bad_x = verify(cigar,ref_seq,qry_seq)
    ref_len,qry_len = len(ref_seq),len(qry_seq)

    print(f"inverted regions flipped: {stats['inverted']} ({stats['inverted_bp']} bp)",file=sys.stderr)
    print(f"blocks chained: {stats['used']} / {stats['blocks']} "
          f"(dropped as fully overlapping: {stats['skipped']})",file=sys.stderr)
    print(f"  alignment from PAF: {stats['from_paf']}, "
          f"reconstructed from SyRI variants: {stats['from_syri']}",file=sys.stderr)
    print(f"  '=' disagreement between PAF and SyRI: {stats['eq_divergence']}",file=sys.stderr)
    print(f"variants dropped as overlapping: {stats['dropped']}",file=sys.stderr)
    print(f"internal skew reconciled: {stats['skew']}",file=sys.stderr)
    print(f"runs realigned by repair pass: {stats['repaired']}",file=sys.stderr)
    print(f"operations: {len(cigar)}",file=sys.stderr)
    for op in "=XID":
        print(f"  {op}: {counts.get(op,0)}",file=sys.stderr)
    print(f"reference consumed: {r} / {ref_len} {'OK' if r == ref_len else 'MISMATCH'}",file=sys.stderr)
    print(f"query consumed:     {q} / {qry_len} {'OK' if q == qry_len else 'MISMATCH'}",file=sys.stderr)
    print(f"'=' positions that are not identical: {bad_eq}",file=sys.stderr)
    print(f"'X' positions that are identical:     {bad_x}",file=sys.stderr)
    print(f"identity: {counts.get('=',0)} / {min(ref_len,qry_len)} "
          f"({100*counts.get('=',0)/min(ref_len,qry_len):.2f}%)",file=sys.stderr)

    text = "".join(f"{length}{op}" for op,length in cigar)
    if out_fname:
        with open(out_fname,'w') as f:
            f.write(text + "\n")
        # the CIGAR aligns against the query with its inversions flipped
        with open(out_fname + ".qry.fa",'w') as f:
            f.write(">inverted_" + read_fasta(qry_fname)[0] + "\n")
            for i in range(0,len(qry_seq),60):
                f.write(qry_seq[i:i+60] + "\n")
    else:
        print(text)
