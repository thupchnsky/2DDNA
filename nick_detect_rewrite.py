#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for nick identification for the rewriting experiment on "GRAINGER". The sequencing files are too large
    to be included here. The link to those files will be provided if requested. There are two parameters eps_1 and eps_2
    that control the number of prefix-suffix pairs we found. Different combination leads to different results. The
    combination here should reproduce the results shown in paper. Note that the detection scheme is slightly different
    from writing and erasure one because we specially design how we assign enzymes into each pool, so we are using a
    different rule suitable for the setting here.

    Usage:
        python nick_detect_rewrite.py
"""

import time
from Bio import SeqIO


# Find the positions of all substrings
def find_all(sub, s):
    index_list = []
    index = s.find(sub)
    while index != -1:
        index_list.append(index)
        index = s.find(sub, index + 1)
    if len(index_list) > 0:
        return index_list
    else:
        return -1


# Count the mismatched positions
def mismatch_count(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


if __name__ == '__main__':
    start = time.time()
    for pool_num in [1, 2, 3, 4, 5, 6, 7, 8]:
        # Sequencing results for rewriting experiment
        if pool_num == 1:
            filePath = '../sequencing_results/rewrite_pool/R1_CGATGTAT_L001_R1_001.fastq'
        elif pool_num == 2:
            filePath = '../sequencing_results/rewrite_pool/R2_TGACCAAT_L001_R1_001.fastq'
        elif pool_num == 3:
            filePath = '../sequencing_results/rewrite_pool/G3_ACAGTGAT_L001_R1_001.fastq'
        elif pool_num == 4:
            filePath = '../sequencing_results/rewrite_pool/N4_GCCAATAT_L001_R1_001.fastq'
        elif pool_num == 5:
            filePath = '../sequencing_results/rewrite_pool/G5_CAGATCAT_L001_R1_001.fastq'
        elif pool_num == 6:
            filePath = '../sequencing_results/rewrite_pool/I6_CTTGTAAT_L001_R1_001.fastq'
        elif pool_num == 7:
            filePath = '../sequencing_results/rewrite_pool/E7_AGTCAACA_L001_R1_001.fastq'
        elif pool_num == 8:
            filePath = '../sequencing_results/rewrite_pool/A8_AGTTCCGT_L001_R1_001.fastq'
        # This is the path to groundtruth oligos in that pool
        poolPath = '../oligo_pools/Level ' + str(pool_num - 1) + '.txt'
        with open(poolPath, 'r') as f:
            ori_lines = f.read().splitlines()

        count = [0] * 6
        nicked_lines = []
        line_count = 1
        for rec in SeqIO.parse(filePath, 'fastq'):
            nicked_lines.append(rec.seq._data)
        # predefined constraints
        eps_1 = 25
        eps_2 = 20
        for line in nicked_lines:
            break_flag = False
            # enzyme 1 Nb.BbvCI
            if line.startswith('TGAGG'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('GCTGAGG', gt_line)
                    if pattern_pos != -1:
                        sub_pos = [i + 2 for i in pattern_pos]
                        for pos in sub_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= eps_2 \
                                            and len(tmp_line) <= pos:
                                        count[0] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 2 Nb.BsmI
            elif line.startswith('CATTC'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('GCATTC', gt_line)
                    if pattern_pos != -1:
                        sub_pos = [i + 1 for i in pattern_pos]
                        for pos in sub_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= eps_2 \
                                            and len(tmp_line) <= pos:
                                        count[1] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 3 Nt.BbvCI
            elif line.startswith('TCAGC'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('CCTCAGC', gt_line)
                    if pattern_pos != -1:
                        sub_pos = [i + 2 for i in pattern_pos]
                        for pos in sub_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1+10 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= eps_2+20 \
                                            and len(tmp_line) <= pos:
                                        count[2] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 4 Nb.BssSI
            elif line.startswith('TCGTG'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('CTCGTG', gt_line)
                    if pattern_pos != -1:
                        sub_pos = [i + 1 for i in pattern_pos]
                        for pos in sub_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= eps_2 \
                                            and len(tmp_line) <= pos:
                                        count[3] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 5 Nt.BspQI
            elif line[-8:-1] == 'GCTCTTC':
                for gt_line in ori_lines:
                    pattern_pos = find_all('GCTCTTC', gt_line)
                    if pattern_pos != -1:
                        for pos in pattern_pos:
                            if mismatch_count(line, gt_line[max(0, pos + 8 - len(line)): pos + 8]) <= eps_1 \
                                    and len(line) <= pos + 8:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[pos + 8: min(pos + 8 + len(tmp_line),
                                                                                     len(gt_line))]) <= eps_2 \
                                            and len(tmp_line) <= len(gt_line) - pos - 8:
                                        count[4] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 6 Nb.BtsI
            elif line.startswith('CACTGC'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('CACTGC', gt_line)
                    if pattern_pos != -1:
                        for pos in pattern_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= eps_2 \
                                            and len(tmp_line) <= pos:
                                        count[5] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
        print(pool_num, ':', count)
        print('Time used:', time.time() - start)
