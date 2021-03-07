#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for Huffman encoding and decoding.
    Usage:
        dictionary, dictionary_inverse = get_dic(input_list)
        s = enc(input_list, dictionary)
        decoded_list = dec_new(s, dictionary_inverse)
"""

from heapq import heappush, heappop, heapify
import numpy as np
import collections 


def enc(l, d):
    # input: list, dictionary
    return ''.join([d[x] for x in l])  # concatenate the binary strings to form a long string


def dec(s, d_i):
    # Huffman decoding when there is no error
    # input: binary string, inverse of the dictionary
    l = []
    while s:
        flag = 0
        for k in d_i:
            if s.startswith(k):
                l.append(d_i[k])
                s = s[len(k):]
                flag = 1
                break
        if flag == 0:
            print('Not found, error!')
    return l


def dec_new(s, d_i):
    # Huffman decoding when there are errors
    # input: binary string, inverse of the dictionary
    l = []
    while s:
        flag = 0
        for k in d_i:
            if s.startswith(k):
                l.append(d_i[k])
                s = s[len(k):]
                flag = 1
                break
        if flag == 1:
            # -2 indicates the end of the information
            if d_i[k] == -2:
                l.pop()
                break
        else:
            # move on to the next bit
            s = s[1:]
    return l


def get_dic(l):
    # generate Huffman dictionary and its inverse based on input
    # input: list
    symb2freq = collections.Counter(l)
    """Huffman encode the given dict mapping symbols to weights"""
    heap = [[wt, [sym, ""]] for sym, wt in symb2freq.items()]
    heapify(heap)
    while len(heap) > 1:
        lo = heappop(heap)
        hi = heappop(heap)
        for pair in lo[1:]:
            pair[1] = '0' + pair[1]
        for pair in hi[1:]:
            pair[1] = '1' + pair[1]
        heappush(heap, [lo[0] + hi[0]] + lo[1:] + hi[1:])
    dictionary = sorted(heappop(heap)[1:], key=lambda p: (len(p[-1]), p))
    data = {x[0]: x[1] for x in dictionary}
    data_inverse = {v: k for k, v in data.items()}
    return data, data_inverse


if __name__ == '__main__':
    l = [1, 1, 2, 3, 4]
    dictionary, dictionary_inverse = get_dic(l)
    s = enc(l, dictionary)
    print('encoded data:', s)
    l_c = dec_new(s, dictionary_inverse)
    print('decode data:', l_c)
    print('original data:', l)

