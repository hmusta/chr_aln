#!/usr/bin/env python3

import sys

OUTDIR=sys.argv[1]
ASSEM_NAME="HG002"
REF_SUFFIX="_PATERNAL"
QRY_SUFFIX="_MATERNAL"

#CHR_RANGE=range(1,23)
CHR_RANGE=range(22,23)
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
eq_total = 0
match_total = 0
for CHR in CHR_RANGE:
    ref_assem, ref_n = fai_get_len(CHR, REF_SUFFIX, CHECK_N)
    qry_assem, qry_n = fai_get_len(CHR, QRY_SUFFIX, CHECK_N)
    minlen_assem = min(ref_assem, qry_assem)

    OUTFILE=f"{OUTDIR}/chr{CHR}_checker.out"
    eq, match, nref, nqry, nlongindelref, nlongindelqry = [int(a) for a in tail(OUTFILE, 2)[0].rstrip().split()]
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

    maxlen = max(ref, qry)

    print(f"{CHR}\t{ref_assem}->{ref}\t{qry_assem}->{qry}\t{eq}\t{100*eq/maxlen:.2f}\t{100*eq/minlen:.2f}")

    ref_total += ref
    qry_total += qry
    eq_total += eq
    match_total += match

min_total = min(ref_total, qry_total)
max_total = max(ref_total, qry_total)
print(f"{ASSEM_NAME}\t{ref_assem_len}->{ref_total}\t{qry_assem_len}->{qry_total}\t{eq_total}\t{100*eq_total/max_total:.2f}\t{100*eq_total/min_total:.2f}")
