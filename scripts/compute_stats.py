#!/usr/bin/env python3

import sys

OUTDIR=sys.argv[1]
ASSEM_NAME="HG002"
REF_SUFFIX="_PATERNAL"
QRY_SUFFIX="_MATERNAL"

CHR_RANGE=range(1,23)
CHECK_N=False

def tail(path, n=2, chunk=4096):
    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        data = b""
        while size > 0 and data.count(b"\n") <= n:
            step = min(chunk, size)
            size -= step
            f.seek(size)
            data = f.read(step) + data
    return data.decode().splitlines()[-n:]

def fai_get_len(CHR, FASUFFIX, CHECK_N=False, FADIR=".", FAEXT="fa"):
    faname=f"chr{CHR}{FASUFFIX}"
    fapath=f"{FADIR}/{faname}.{FAEXT}"
    path=fapath + ".fai"
    tot_len = 0
    with open(path, 'r') as f:
        line = f.readlines()
        assert(len(line) == 1)
        line = line[0].rstrip().split("\t")
        assert(line[0] == faname)
        tot_len = int(line[1])

    tot_n = 0
    if CHECK_N:
        tot_len_fa = 0
        with open(fapath, 'r') as f:
            for line in f:
                if line[0] == ">":
                    continue
                tot_len_fa += len(line) - 1
                tot_n += line.count('N')
            assert(tot_len_fa == tot_len)
    return tot_len, tot_n

ref_assem_len = 0
qry_assem_len = 0
ref_total = 0
qry_total = 0
ref_longindel_total = 0
qry_longindel_total = 0
eq_total = 0
match_total = 0
for CHR in CHR_RANGE:
    ref_assem, ref_n = fai_get_len(CHR, REF_SUFFIX, CHECK_N)
    qry_assem, qry_n = fai_get_len(CHR, QRY_SUFFIX, CHECK_N)
    minlen_assem = min(ref_assem, qry_assem)

    OUTFILE=f"{OUTDIR}/chr{CHR}_checker.out"
    try:
        eq, match, nref, nqry, nlongindelref, nlongindelqry = [int(a) for a in tail(OUTFILE, 2)[0].rstrip().split()]
    except ValueError:
        print(f"Skipping chr{CHR} since it's incomplete",file=sys.stderr)
        continue
    if CHECK_N and (nref != ref_n or nqry != qry_n):
        print(f"chr{CHR}:\t{nref} != {ref_n} or {nqry} != {qry_n}")
        print(eq,match,nref,nqry,nlongindelref,nlongindelqry)
    assert(not CHECK_N or nref == ref_n)
    assert(not CHECK_N or nqry == qry_n)
    assert(eq <= minlen_assem)
    assert(match <= minlen_assem)
    assert(nref + nlongindelref <= ref_assem)
    assert(nqry + nlongindelqry <= qry_assem)

    ref_assem_len += ref_assem
    qry_assem_len += qry_assem

    ref = ref_assem - nref
    qry = qry_assem - nqry
    minlen = min(ref, qry)
    assert(eq <= minlen)
    assert(match <= minlen)
    assert(nlongindelref <= ref)
    assert(nlongindelqry <= qry)

    ref_longindel_total += nlongindelref
    qry_longindel_total += nlongindelqry

    maxlen = max(ref, qry)

    eq_pc_min = 100*eq/maxlen
    eq_pc_max = 100*eq/minlen
    match_pc_min = 100*match/maxlen
    match_pc_max = 100*match/minlen

    minlen_nolong = min(ref-nlongindelref,qry-nlongindelqry)
    maxlen_nolong = max(ref-nlongindelref,qry-nlongindelqry)
    eq_nolong_min = 100*eq/maxlen_nolong
    eq_nolong_max = 100*eq/minlen_nolong

    eq_pc = 100*eq/match

    print(f"{CHR}\t{ref_assem}->{ref}\t{qry_assem}->{qry}\t{eq}\t{eq_pc_min:.2f}-{eq_pc_max:.2f}\t{match_pc_min:.2f}-{match_pc_max:.2f}\t{eq_pc:.2f}\t{eq_nolong_min:.2f}-{eq_nolong_max:.2f}")

    ref_total += ref
    qry_total += qry
    eq_total += eq
    match_total += match

minlen = min(ref_total, qry_total)
maxlen = max(ref_total, qry_total)
eq_pc_min = 100*eq_total/maxlen
eq_pc_max = 100*eq_total/minlen
match_pc_min = 100*match_total/maxlen
match_pc_max = 100*match_total/minlen

minlen_nolong = min(ref_total-ref_longindel_total,qry_total-qry_longindel_total)
maxlen_nolong = max(ref_total-ref_longindel_total,qry_total-qry_longindel_total)
eq_nolong_min = 100*eq_total/maxlen_nolong
eq_nolong_max = 100*eq_total/minlen_nolong

eq_pc = 100*eq_total/match_total

print(f"{ASSEM_NAME}\t{ref_assem_len}->{ref_total}\t{qry_assem_len}->{qry_total}\t{eq_total}\t{eq_pc_min:.2f}-{eq_pc_max:.2f}\t{match_pc_min:.2f}-{match_pc_max:.2f}\t{eq_pc:.2f}\t{eq_nolong_min:.2f}-{eq_nolong_max:.2f}")
