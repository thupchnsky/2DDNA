#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for all conversion steps needed in proposed encoding-decoding scheme.
"""

import math
import utils.hilbert as hilbert
import numpy as np
import os


def phi(b1, b2):
    # convert two binary string into one DNA sequence
    d = list('ATCG')
    b = list(np.binary_repr(i, width=2) for i in range(4))
    b2d = {k: v for k, v in zip(b, d)}
    return ''.join([b2d[x+y] for x, y in zip(b1, b2)])  # '00':'A','01':'T','10':'C','11':'G'


def phi_r(dna):
    # reverse mapping
    d = list('ATCG')
    b = list(np.binary_repr(i, width=2) for i in range(4))
    d2b = {k: v for k, v in zip(d, b)}
    s = [d2b[x] for x in dna]
    b_1 = [x[0] for x in s]
    b_2 = [x[1] for x in s]
    return [''.join(b_1), ''.join(b_2)]


def DNA_to_level(read):
    # Here read is an Bio record
    # read is given format ['file', '0', 'color', '0', 'level', '0', 'block', '3']
    name = read.id
    name = name.split('_')
    return int(name[5])


def quantize(l, q):
    # input: l is a list and q is the number of bits to represent each value
    n = 2 ** q
    return [math.floor(x*n/256) for x in l]  # quantized values


def dequantize(l, q):
    # input: l is a list and q is the number of bits to represent each value
    n = 2 ** q
    return [int(x * 256/n) + int(128/n) for x in l]


def dequantize_new(l, q):
    # Instead of using mid value, use a random offset
    # input: l is a list and q is the number of bits to represent each value
    n = 2 ** q
    return [int(x * 256 / n) + np.random.randint(32) for x in l]


def sub_dequantize(l, q, levels):
    # only dequantize specific levels, mask others as white
    # input: l is a list and q is the number of bits to represent each value
    n = 2 ** q
    T = []
    for x in l:
        if x in levels:
            T.append(int(x * 256/n) + int(128/n))
        else:
            T.append(255)
    return T


def matrix_2_vector(m, path):
    # use hilbert curve to flatten the matrix
    x = len(m[0])  # number of columns
    y = len(m)  # number of rows
    filename = path + 'hilbert_curve_' + str(y) + '_' + str(x) + '_dictionary.npy'
    if not os.path.exists(filename):
        h = hilbert.hilbert_get(x, y)
        np.save(filename, h)
    else:
        h = np.load(filename)
    return [m[i[0]][i[1]] for i in h]  # flatten the matrix: m[0][0], m[1][0],...


def vector_2_level(c, f):
    # c is the color, f is a vector
    unq = list(set(f))
    unq.sort()
    l = []
    for i in unq:
        T = [j for j, x in enumerate(f) if x == i]
        l.append([c, i, T])
    return l  # [1, 1, 2 , 1,..] ->{1: [0,1,3,...], 2:[...]}


def level_2_vector(l):
    # l is a list of different level position
    c = l[0][0]
    total_len = 0
    for i in range(len(l)):
        total_len += (len(l[i]) - 2)
    f = [1] * total_len
    for i in range(len(l)):
        for j in range(2,len(l[i])):
            f[l[i][j]] = l[i][1]
    return c, f


def level_2_vector_new(l):
    total_len = 0
    for item in l:
        total_len += len(item)
    f = [1] * total_len
    for idx in range(len(l)):
        for item in l[idx]:
            f[min(total_len - 1, item)] = idx
    return f


def diff_enc(l):
    # Differential encoding. Add -1 as synchronizing markers
    result_l = [-1, l[0]]
    for i in range(1, len(l)):
        if i % 30 == 0:
            result_l += [-1, l[i]]
        else:
            result_l += [l[i] - l[i-1]]
    return result_l


def diff_dec(l):
    # Differential decoding. Update value when -1 appears
    result_l = [l[1]]
    i = 2
    while i < len(l):
        if l[i] == -1:
            if i + 1 < len(l):
                result_l += [l[i+1]]
                i = i + 2
            else:
                return result_l
        else:
            result_l += [l[i] + result_l[-1]]
            i += 1
    return result_l


def decide_level(dna_seq, pre_dict, suf_dict):
    pre_len = len(pre_dict['0'])
    suf_len = len(suf_dict['0'])
    pre_seq = dna_seq[0:pre_len]
    suf_seq = dna_seq[-suf_len:]
    err_bits = []
    for i in range(len(pre_dict)):
        err_bits.append(str_cmp(pre_seq, pre_dict[str(i)]) + str_cmp(suf_seq, suf_dict[str(i)]))
    return err_bits.index(min(err_bits))


def str_cmp(s1, s2):
    count = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count += 1
    return count


if __name__ == '__main__':
    for _ in range(1000):
        l = list(np.random.randint(100, size=1000))
        l.sort()
        a = diff_enc(l)
        b = diff_dec(a)
        if l != b:
            print('Error happens in differential encoding!')
            print(l)
            print(a)
            print(b)

