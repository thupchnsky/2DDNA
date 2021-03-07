#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for nick identification for the erasure experiment on "ILLINOIS". The sequencing files are too large
    to be included here. The link to those files will be provided if requested. The decision rule here should be the
    same as first round writing experiments for fair comparison.

    Usage:
        python nick_detect_erase.py
"""

import time
from Bio import SeqIO


# find the positions of all substrings
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


# count the mismatched positions
def mismatch_count(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


if __name__ == '__main__':
    start = time.time()
    for pool_num in [1, 2, 3, 4, 5, 6, 7, 8]:
        # Sequencing results for erasure experiment
        if pool_num == 1:
            filePath = '../sequencing_results/erase_pool/I1_L_CAGATCAT_L001_R1_001.fastq'
        elif pool_num == 2:
            filePath = '../sequencing_results/erase_pool/L2_L_CTTGTAAT_L001_R1_001.fastq'
        elif pool_num == 3:
            filePath = '../sequencing_results/erase_pool/L3_L_AGTCAACA_L001_R1_001.fastq'
        elif pool_num == 4:
            filePath = '../sequencing_results/erase_pool/I4_L_AGTTCCGT_L001_R1_001.fastq'
        elif pool_num == 5:
            filePath = '../sequencing_results/erase_pool/N5_L_ATGTCAGA_L001_R1_001.fastq'
        elif pool_num == 6:
            filePath = '../sequencing_results/erase_pool/O6_L_CCGTCCCG_L001_R1_001.fastq'
        elif pool_num == 7:
            filePath = '../sequencing_results/erase_pool/I7_L_GTCCGCAC_L001_R1_001.fastq'
        elif pool_num == 8:
            filePath = '../sequencing_results/erase_pool/S8_L_GTGAAACG_L001_R1_001.fastq'
        # This is the path to groundtruth oligos in that pool
        poolPath = '../oligo_pools/Level ' + str(pool_num - 1) + '.txt'
        with open(poolPath, 'r') as f:
            ori_lines = f.read().splitlines()

        count = [0] * 7
        nicked_lines = []
        line_count = 1
        for rec in SeqIO.parse(filePath, 'fastq'):
            nicked_lines.append(rec.seq._data)

        eps_1 = 9
        eps_2 = 9
        for line in nicked_lines:
            break_flag = False
            # enzyme 1 Nb.BtsI
            if line.startswith('CACTGC'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('CACTGC', gt_line)
                    if pattern_pos != -1:
                        for pos in pattern_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= len(tmp_line)/3 \
                                            and len(tmp_line) <= pos:
                                        count[0] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 2 Nt.BstNBI
            elif line.endswith('GAGTC'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('GAGTC', gt_line)
                    if pattern_pos != -1:
                        for pos in pattern_pos:
                            if mismatch_count(line, gt_line[max(0, pos + 5 - len(line)): pos + 5]) <= eps_1 \
                                    and len(line) <= pos + 5:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[pos + 5: min(pos + 5 + len(tmp_line), len(gt_line))]) <= len(tmp_line)/3 \
                                            and len(tmp_line) <= len(gt_line) - pos - 5:
                                        count[1] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 3 Nb.BssSI
            elif line.startswith('TCGTG'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('CTCGTG', gt_line)
                    if pattern_pos != -1:
                        sub_pos = [i + 1 for i in pattern_pos]
                        for pos in sub_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= len(tmp_line)/3 \
                                            and len(tmp_line) <= pos:
                                        count[2] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 4 Nt.AIwI
            elif line[-9:-4] == 'GGATC':
                for gt_line in ori_lines:
                    pattern_pos = find_all('GGATC', gt_line)
                    if pattern_pos != -1:
                        for pos in pattern_pos:
                            if mismatch_count(line, gt_line[max(0, pos + 9 - len(line)): pos + 9]) <= eps_1 \
                                    and len(line) <= pos + 9:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[pos + 9: min(pos + 9 + len(tmp_line), len(gt_line))]) <= len(tmp_line)/3 \
                                            and len(tmp_line) <= len(gt_line) - pos - 9:
                                        count[3] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 5 Nt.BsmAI
            elif line[-6:-1] == 'GTCTC':
                for gt_line in ori_lines:
                    pattern_pos = find_all('GTCTC', gt_line)
                    if pattern_pos != -1:
                        for pos in pattern_pos:
                            if mismatch_count(line, gt_line[max(0, pos + 6 - len(line)): pos + 6]) <= eps_1 \
                                    and len(line) <= pos + 6:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[pos + 6: min(pos + 6 + len(tmp_line), len(gt_line))]) <= len(tmp_line)/2.5 \
                                            and len(tmp_line) <= len(gt_line) - pos - 6:
                                        count[4] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 6 Nb.BsmI
            elif line.startswith('CATTC'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('GCATTC', gt_line)
                    if pattern_pos != -1:
                        sub_pos = [i + 1 for i in pattern_pos]
                        for pos in sub_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= len(tmp_line)/3 \
                                            and len(tmp_line) <= pos:
                                        count[5] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
            # enzyme 7 Nb.BsrDI
            elif line.startswith('CATTGC'):
                for gt_line in ori_lines:
                    pattern_pos = find_all('CATTGC', gt_line)
                    if pattern_pos != -1:
                        for pos in pattern_pos:
                            if mismatch_count(line, gt_line[pos: min(pos + len(line), len(gt_line))]) <= eps_1 \
                                    and len(line) <= len(gt_line) - pos:
                                for tmp_line in nicked_lines:
                                    if mismatch_count(tmp_line, gt_line[max(0, pos - len(tmp_line)): pos]) <= len(tmp_line)/3 \
                                            and len(tmp_line) <= pos:
                                        count[6] += 1
                                        break_flag = True
                                        break
                                break
                    if break_flag:
                        break
        print(pool_num, ':', count)
        print('Time used:', time.time() - start)
