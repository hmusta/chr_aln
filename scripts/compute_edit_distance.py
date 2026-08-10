#!/usr/bin/env python3

import sys
import itertools
import pathlib
import numpy as np

REF_CONSUME="D"
QRY_CONSUME="I"
edits = ["X",REF_CONSUME,QRY_CONSUME,"S","H"]
consume_ref = ["M","=","X",REF_CONSUME]
consume_query = ["M","=","X",QRY_CONSUME]

_COMPLEMENT = str.maketrans(
    "ACGTNacgtnRYSWKMBDHVryswkmbdhv",
    "TGCANtgcanYRSWMKVHDByrwsmkvhdb",
)

indel_length_cutoff = 25

total_length_a = 0
total_length_b = 0
total_full_length_a = 0
total_full_length_b = 0
total_dist = 0
total_id = 0
total_matches = 0
total_length_short_indels_a = 0
total_length_short_indels_b = 0
nchrom = 0

#print("chr\tdist\tid\tlength a\tlength b\tid frac. match\tid frac. sh.\tid frac.\tmax diag.")
#print("chr\tid frac. match\tid frac. sh.\tid frac.\tmax diag.\trefL\trefN\tqryL\tqryN")

check_N = bool(int(sys.argv[1]))
check_eq = bool(int(sys.argv[2]))
print(check_eq,check_N,file=sys.stderr)
nins = 4
args_iter = list(itertools.zip_longest(*[iter(sys.argv[3:])]*nins))
for ref_fname, query_fname, fname, inv_bed in args_iter:
    name = fname.split("/")[-1].split(".")[0]
    print(name,file=sys.stderr)
    with open(ref_fname, 'r') as f:
        ref_header = f.readline().rstrip().split()[0][1:]
        ref_seq = "".join(line.rstrip().upper() for line in f.readlines())
    with open(query_fname, 'r') as f:
        query_header = f.readline().rstrip().split()[0][1:]
        query_seq = "".join(line.rstrip().upper() for line in f.readlines())
    ref_inv = np.zeros(len(ref_seq), dtype=bool)
    qry_inv = np.zeros(len(query_seq), dtype=bool)
    if pathlib.Path(inv_bed).resolve().is_file():
        with open(inv_bed, 'r') as f:
            for line in f:
                line = line.rstrip().split("\t")
                assert(line[0] == ref_header)
                assert(line[6] == query_header)
                ref_begin = int(line[1])
                ref_end = int(line[2])
                ref_inv[ref_begin:ref_end] = True
                query_begin = int(line[7])
                query_end = int(line[8])
                qry_inv[query_begin:query_end] = True
                qf = query_seq[:query_begin]
                qi = query_seq[query_begin:query_end].translate(_COMPLEMENT)[::-1]
                ql = query_seq[query_end:]
    with open(fname, 'r') as f:
        line = f.readline()
        if line == '':
            continue
        cigar = line.rstrip().split()[-1]
        cur_num = 0
        length = 0
        length_q = 0
        short_indels_a = 0
        short_indels_b = 0
        dist = 0
        cur_id = 0
        cur_matches = 0
        cur_diag = 0
        max_diag = 0
        ref_N = 0
        query_N = 0
        for c in cigar:
            if c.isdigit():
                cur_num = cur_num * 10 + int(c)
            else:
                if c in edits:
                    dist += cur_num
                if c == QRY_CONSUME:
                    cur_diag -= cur_num
                    cur_N = query_seq[length_q:length_q+cur_num].count('N')
                    if abs(cur_diag) > abs(max_diag):
                        max_diag = cur_diag
                    if cur_num-cur_N <= indel_length_cutoff:
                        short_indels_b += cur_num-cur_N
                    query_N += cur_N
                    length_q += cur_num
                elif c == REF_CONSUME:
                    cur_diag += cur_num
                    cur_N = ref_seq[length:length+cur_num].count('N')
                    if abs(cur_diag) > abs(max_diag):
                        max_diag = cur_diag
                    if cur_num-cur_N <= indel_length_cutoff:
                        short_indels_a += cur_num-cur_N
                    ref_N += cur_N
                    length += cur_num
                elif c == "M":
                    assert(not check_N)
                    assert(np.all(ref_inv[length:length+cur_num] == qry_inv[length_q:length_q+cur_num]))
                    ref_w = ref_seq[length:length+cur_num]
                    qry_w = query_seq[length_q:length_q+cur_num]
                    num_N = sum(rr == 'N' or qq == 'N' for rr,qq in zip(ref_w,qry_w))
                    assert(cur_num == sum(rr == 'N' or qq == 'N' for rr,qq in zip(ref_w,qry_w)))
                    ref_N += ref_w.count('N')
                    query_N += qry_w.count('N')
                elif c == "=":
                    assert(np.all(ref_inv[length:length+cur_num] == qry_inv[length_q:length_q+cur_num]))
                    ref_w = ref_seq[length:length+cur_num]
                    qry_w = query_seq[length_q:length_q+cur_num]
                    num_N = sum(rr == 'N' or qq == 'N' for rr,qq in zip(ref_w,qry_w))
                    assert(num_N <= cur_num)
                    num_obs_eq = sum(rr != 'N' and qq != 'N' and rr == qq for rr,qq in zip(ref_w,qry_w))
                    num_obs_neq = sum(rr != 'N' and qq != 'N' and rr != qq for rr,qq in zip(ref_w,qry_w))
                    if num_obs_neq > 0:
                        print(f"{length+1}-{length+cur_num}",file=sys.stderr)
                        print(f"{length_q+1}-{length_q+cur_num}",file=sys.stderr)
                        print(c,"\t","=",num_obs_eq,"X",num_obs_neq,"N",num_N,"total",cur_num,file=sys.stderr)
                        print(ref_w,file=sys.stderr)
                        print(qry_w,file=sys.stderr)
                        print(file=sys.stderr)
                    assert(not check_eq or num_obs_neq == 0)

                    num_eq = num_obs_eq
                    cur_id += num_eq
                    cur_matches += num_eq
                    ref_N += ref_w.count('N')
                    query_N += qry_w.count('N')
                    length += cur_num
                    length_q += cur_num
                elif c == "X":
                    assert(np.all(ref_inv[length:length+cur_num] == qry_inv[length_q:length_q+cur_num]))
                    ref_w = ref_seq[length:length+cur_num]
                    qry_w = query_seq[length_q:length_q+cur_num]

                    num_N = sum(rr == 'N' or qq == 'N' for rr,qq in zip(ref_w,qry_w))
                    assert(check_N or num_N <= cur_num)

                    num_obs_neq = sum(rr != 'N' and qq != 'N' and rr != qq for rr,qq in zip(ref_w,qry_w))
                    num_obs_eq = sum(rr != 'N' and qq != 'N' and rr == qq for rr,qq in zip(ref_w,qry_w))
                    assert(not check_eq or num_obs_eq == 0)

                    num_neq = num_obs_neq - num_N
                    cur_matches += num_neq
                    ref_N += ref_w.count('N')
                    query_N += qry_w.count('N')
                    length += cur_num
                    length_q += cur_num
                assert(length <= len(ref_seq))
                assert(length_q <= len(query_seq))
                assert(cur_diag == length - length_q)
                cur_num = 0
    assert(length == len(ref_seq))
    assert(length_q == len(query_seq))
    final_diag = length - length_q
    if length_q > length:
        final_diag *= -1
        max_diag *= -1
    full_length = length
    full_length_q = length_q
    length -= ref_N
    length_q -= query_N
    id_sh_indels_a = cur_id / (cur_matches + short_indels_a)
    id_sh_indels_b = cur_id / (cur_matches + short_indels_b)
    id_sh_indels_min = min(id_sh_indels_a,id_sh_indels_b)
    id_sh_indels_max = max(id_sh_indels_a,id_sh_indels_b)

    id_frac_a = cur_id / length
    id_frac_b = cur_id / length_q
    id_frac_min = min(id_frac_a,id_frac_b)
    id_frac_max = max(id_frac_a,id_frac_b)
    minlen = min(length,length_q)
    maxlen = max(length,length_q)
    print(f"{name[3:]}\t{full_length}->{length}\t{full_length_q}->{length_q}\t{cur_id}\t{100*cur_id/maxlen:.2f}\t{100*cur_id/minlen:.2f}")
    #print(f"{name[3:].rjust(3)}\t{(cur_id / cur_matches):.5f}\t{id_sh_indels_min:.5f}-{id_sh_indels_max:.5f}\t{id_frac_min:.5f}-{id_frac_max:.5f}\t{max_diag} / {final_diag}\t{length}\t{ref_N}\t{length_q}\t{query_N}")
    #print(f"{name[3:].rjust(3)}\t{dist}\t{cur_id}\t{length}\t{length_q}\t{(cur_id / cur_matches):.3f}\t{id_sh_indels_min:.3f}-{id_sh_indels_max:.3f}\t{id_frac_min:.3f}-{id_frac_max:.3f}\t{max_diag} / {final_diag}")

    total_full_length_a += full_length
    total_full_length_b += full_length_q
    total_dist += dist
    total_length_a += length
    total_length_b += length_q
    total_id += cur_id
    total_matches += cur_matches
    total_length_short_indels_a += cur_matches + short_indels_a
    total_length_short_indels_b += cur_matches + short_indels_b
    nchrom += 1

if nchrom > 1:
    fullname="HG002"
    minlen=min(total_length_a,total_length_b)
    maxlen=max(total_length_a,total_length_b)
    print(f"{fullname}\t{total_full_length_a}->{total_length_a}\t{total_full_length_b}->{total_length_b}\t{total_id}\t{100*total_id/maxlen:.2f}\t{100*total_id/minlen:.2f}")
#    print(
#        f"total dist: {total_dist}\n" + \
#        f"total id: {total_id}\n" + \
#        f"total matches: {total_matches}\n" + \
#        f"total length a: {total_length_a}\n" + \
#        f"total length b: {total_length_b}\n" + \
#        f"id frac. match: {(total_id / total_matches):.5f}\n" + \
#        f"id frac. sh. a: {(total_id / total_length_short_indels_a):.5f}\n" + \
#        f"id frac. sh. b: {(total_id / total_length_short_indels_b):.5f}\n" + \
#        f"id frac. a: {(total_id / total_length_a):.5f}\n" + \
#        f"id frac. b: {(total_id / total_length_b):.5f}" \
#        )
