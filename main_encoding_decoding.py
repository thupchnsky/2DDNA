#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for image encoding and decoding. Configuration stored in file "input.json" includes:
    prefix: prefix used for each intensity level pool
    suffix: suffix used for each intensity level pool
    channel: 3 nts DNA block used for each color channel (B, G, R)
    bit:
        file_id: number of bits used to store file index
        block_id: number of bits used to store block index
        q: number of bits used to represent intensity level after quantization
        block: [number of bits in one information sub-block, number of information sub-blocks in one oligo]
    gc:
        id: [length of DNA block used to represent binary address block, min number of symbol G/C in this DNA block,
        max number of symbol G/C in this DNA block]
        block: [length of DNA block used to represent binary information sub-block, min number of symbol G/C in this
        DNA block, max number of symbol G/C in this DNA block]
    path:
        input: list of path for all target images
        output: path where all results are stored
    order:
        levels: which pools we should store the results
    size: image sizes

    Usage:
        if both encoding and decoding is needed, then use
            python main_encoding_decoding.py --config="input.json"
        if only decoding is needed (when consensus sequencing reads are obtained), then use
            python main_encoding_decoding.py --config="input.json" --decoding_only=True --read_path=PATH_TO_READS
        One example of PATH_TO_READS is included in "reads_example/consensus_seq_with_errors.txt"
        Note that one round of encoding is needed before any kind of decoding, since the Hilbert curve and Huffman
        dictionary needs to be generated.
"""

# import built-in packages
import numpy as np
import json
import cv2
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pickle
import argparse
# import our packages
import utils.huffman as huffman
import utils.conversion as conversion


def str_cmp(s1, s2):
    """
    Compute the Hamming distance between two strings of the same length.
    """
    count = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count += 1
    return count


def json_load(path):
    """
    Load information stored in json file.
    """
    with open(path) as data_file:  # open & read json
        return json.load(data_file)


def json_save(obj, name):
    """
    Save information to json file.
    """
    with open(name, 'w') as fout:
        json.dump(obj, fout)


def decide_color(str1, ref1, ref2, ref3):
    """
    Choose the index of ref sequence that have the smallest Hamming distance with str1
    """
    count = [str_cmp(str1, ref1), str_cmp(str1, ref2), str_cmp(str1, ref3)]
    return count.index(min(count))


class gc(object):
    """
        Generate the constrained mapping between DNA block and binary string. Note that not all (l, gc_min, gc_max)
        triple can generate valid mappings.

        Attributes:
            l: length of the mapped DNA block
            gc_min: min number of symbol G/C in this DNA block
            gc_max: max number of symbol G/C in this DNA block

        Usage:
            gc_blck = gc(10, 4, 6)
    """
    def __init__(self, l, gc_min, gc_max):
        self.l, self.gc_min, self.gc_max = l, gc_min, gc_max
        self.set_c_1()

    def run_1(self, seq):
        # return the max run length of '1' in string
        return max(map(len, seq.split('0')))

    def set_c_1(self):
        # split possible binary strings into two parts of different length. Map shorter part to another binary string of
        # the same length of the longer part to ensure the sparsity of '1'.
        # return the mapping and the reverse mapping
        v = []
        for i in range(pow(2, self.l)):
            bin_seq = np.binary_repr(i, width=self.l)
            if self.gc_min <= bin_seq.count('1') <= self.gc_max and self.run_1(bin_seq) < 3 and bin_seq[0:2] != '11':
                v.append(bin_seq)
        self.l_1 = int(np.log2(len(v)))
        v = v[0: pow(2, self.l_1)]
        k = list(np.binary_repr(i, width=self.l_1) for i in range(pow(2, self.l_1)))

        self.c_1 = {x: y for x, y in zip(k, v)}
        self.c_1_r = {x: y for x, y in zip(v, k)}

    def bit2dna(self, b):
        # map binary string to DNA block
        if len(b) > self.l + self.l_1:
            print('Long string')
            return False
        elif len(b) < self.l + self.l_1:
            print('Short string')
            return False
        else:
            b_1 = self.c_1[b[0:self.l_1]]
            b_2 = b[self.l_1:]
            return conversion.phi(b_1, b_2)

    def dna2bit(self, d):
        # map DNA block back to binary string. Errors can exist.
        [b1, b2] = conversion.phi_r(d)
        if b1 in self.c_1_r.keys():
            # if there is no error
            b = self.c_1_r[b1]
        else:
            # if error exists
            min_key = ''
            min_count = 100
            for key in self.c_1_r:
                count = str_cmp(key, b1)
                if count < min_count:
                    min_count = count
                    min_key = key
            b = self.c_1_r[min_key]
        return b + b2


class channel(object):
    """
        encode binary information into pools of oligos.

        Attributes:
            data: General setting read from configuration file
            gc_blck: gc object for information sub-blocks
            gc_adrs: gc object for address block

        Usage:
            records, reads_dna, level_size = channel(data).enc(binary_data)
    """

    def __init__(self, data):
        # define two mappings. One for address block and one for information sub-blocks
        self.data = data
        self.gc_blck = gc(data['gc']['block'][0], data['gc']['block'][1], data['gc']['block'][2])
        self.gc_adrs = gc(data['gc']['id'][0], data['gc']['id'][1], data['gc']['id'][2])

    def enc(self, bgr):
        # encode binary information bgr into pools of oligos reads_dna
        # input: [file number, color, level, binary string]
        records = []  # initialize DNA sequences
        reads_dna = []
        for x in bgr:  # loop over vectors
            rds, tmp_dna = self.form_block(x)
            records += rds  # add the new reads
            reads_dna += tmp_dna

        SeqIO.write((x for x in records), self.data['path']['output'] + 'reads.fasta', "fasta")  # write the records into file

        # store oligos by levels
        levels = list(set([conversion.DNA_to_level(read) for read in records]))
        levels.sort()
        level_size = []
        for level in levels:
            read_filter = [read for read in records if conversion.DNA_to_level(read) == level]  # filter only that level
            SeqIO.write((x for x in read_filter), self.data['path']['output'] + 'level_' + str(level) + '.fasta', "fasta")  # write the records into file
            level_size.append(len(read_filter))

        return records, reads_dna, level_size  # return generated oligos and related information

    def form_block(self, b):
        # concatenating different DNA blocks to obtain oligo of the same length
        # input: [file number, color, level, binary string]
        # information block
        add = np.random.randint(2, size=(-len(b[3])) % (self.data['bit']['block'][0] * self.data['bit']['block'][1]))  # generate extra random bits to fit blocks
        m_b = b[3] + ''.join(map(str, add))  # modify binary sequence by adding extra bits
        m_b = [m_b[i:i + self.data['bit']['block'][0]] for i in range(0, len(m_b), self.data['bit']['block'][0])]  # m_b divide it into blocks of length self.data['bit']['block'][0]
        dna = [self.gc_blck.bit2dna(b) for b in m_b]  # dna string
        dna = ''.join(dna)  # form a long string
        infs = [dna[i:i + self.data['bit']['block'][1] * self.data['gc']['block'][0]] for i in range(0, len(dna), self.data['bit']['block'][1] * self.data['gc']['block'][0])]  # blocks of length self.data['bit']['block'][1] * self.data['gc']['block'][0]

        # id block
        color = self.data['channel'][b[1]]
        file_id = np.binary_repr(b[0], width=self.data['bit']['file_id'])
        level_id = np.binary_repr(b[2], width=self.data['bit']['q'])
        blks_id = [np.binary_repr(i, width=self.data['bit']['block_id']) for i in range(len(infs))]
        ids = [color + self.gc_adrs.bit2dna(file_id + b + level_id) for b in blks_id]  # c + file_id + block_id + level_id
        dna = [self.data['prefix'][str(b[2])] + x + y + self.data['suffix'][str(b[2])] for x, y in zip(ids, infs)]  # Assembling the blocks

        # records
        records = []
        for i, d in enumerate(dna):
            r = SeqRecord(Seq(d, IUPAC.protein), id='file_' + str(b[0]) + '_color_' + str(b[1]) + '_level_' + str(b[2]) + '_block_' + str(i), name='', description='')
            records.append(r)
        return records, dna


def enc(f, path_in, path_out, q):
    """
    Encode one image into binary strings using proposed encoding scheme.
    :param f: int, file index
    :param path_in: str, path of the image
    :param path_out: str, path of the folder storing results
    :param q: number of bits used for quantization
    :return: bgr: list, encoded binary strings; bgr_q: list, dequantized pixel values which can be used for
             reconstruction
    """
    # load image
    img = cv2.imread(path_in)  # read image
    row, col = img.shape[0:2]  # get image information

    # channels, quantization and scale
    bgr = cv2.split(img)  # blue, green, red channels
    bgr = [[conversion.quantize(x, q) for x in y] for y in bgr]  # quantize and scale matrix

    # save quantized image
    bgr_q = [[conversion.dequantize(x, q) for x in y] for y in bgr]  # quantize matrix
    bgr_q = [np.array(x, dtype=np.uint8) for x in bgr_q]  # convert list to array
    img_q = cv2.merge(bgr_q)  # merge channels
    cv2.imwrite(path_out + 'file_' + str(f) + '_bits_' + str(q) + '.png', img_q)  # save image

    # encoding
    bgr = [[i, conversion.matrix_2_vector(x, path_out)] for i, x in enumerate(bgr)]  # matrix to vector
    bgr = [conversion.vector_2_level(x[0], x[1]) for x in bgr]  # vector to level
    bgr = bgr[0] + bgr[1] + bgr[2]  # [color, level, vector]
    bgr = [[x[0], x[1], conversion.diff_enc(x[2]) + [-2]] for x in bgr]  # differential encoding

    # construct huffman dictionary
    l = []  # initialize list
    for x in bgr:  # loop over vectors
        l = l + x[2]

    dictionary, dictionary_inverse = huffman.get_dic(l)  # get huffman dictionary
    json_save(dictionary, path_out + 'file_' + str(f) + '_huffman_dictionary.json')  # save huffman dictionary
    json_save(dictionary_inverse, path_out + 'file_' + str(f) + '_huffman_dictionary_inverse.json')

    # huffman encoding and add file index
    bgr = [[f, x[0], x[1], huffman.enc(x[2], dictionary)] for x in bgr]  # [file_id, color_id, level_id, binary vector]
    return bgr, bgr_q


def enc_dna(data):
    """
    Encode all images into DNA blocks using proposed encoding scheme.
    :param data: dict, general setting read from configuration file
    :return: records: list, Bio records; reads_dna: list, encoded DNA sequences; ori_q: list, dequantized information;
             level_size: list, number of oligos for each intensity level pool
    """
    # copy the input into the output folder
    json_save(data, data['path']['output'] + '/input.json')

    binary_data = []  # initialize
    ori_q = []
    for f, path in enumerate(data['path']['input']):  # f: file number
        file_binary, file_bgr_q = enc(f, path, data['path']['output'], data['bit']['q'])
        binary_data = binary_data + file_binary
        ori_q += [file_bgr_q]
        print('file ' + str(f) + ' encoding finished.')

    records, reads_dna, level_size = channel(data).enc(binary_data)
    return records, reads_dna, ori_q, level_size


def dna2levels(reads_dna, data):
    """
    Converting DNA oligos into binary strings organized by levels
    :param reads_dna: list, consensus reads of DNA oligos
    :param data: dict, general setting read from configuration file
    :return: binary_seq, list, decoded binary strings
    """
    # initialize parameters
    pre_len = len(data['prefix']['0'])
    suf_len = len(data['suffix']['0'])
    color_len = 3
    adrs_len = data['gc']['id'][0]
    blks_len = data['gc']['block'][0]

    # count number of levels for each color and how many blocks for each level each color for concatenation
    max_level = [-1] * 3
    level_blocks = [[], [], []]
    gc_blck = gc(data['gc']['block'][0], data['gc']['block'][1], data['gc']['block'][2])
    gc_adrs = gc(data['gc']['id'][0], data['gc']['id'][1], data['gc']['id'][2])

    for line in reads_dna:
        adrs = gc_adrs.dna2bit(line[pre_len + color_len: pre_len + color_len + adrs_len])
        level_id = int(adrs[-data['bit']['q']:], 2)
        color_id = decide_color(line[pre_len: pre_len + color_len], data['channel'][0], data['channel'][1],
                                data['channel'][2])
        if level_id > max_level[color_id]:
            level_blocks[color_id] += [0] * (level_id - max_level[color_id])
            level_blocks[color_id][level_id] = 1
            max_level[color_id] = level_id
        else:
            level_blocks[color_id][level_id] += 1

    binary_seq = [[], [], []]
    for i in range(3):
        for j in range(len(level_blocks[i])):
            binary_seq[i] += [[''] * level_blocks[i][j]]

    # put oligos into blocks
    for line in reads_dna:
        adrs = gc_adrs.dna2bit(line[pre_len + color_len: pre_len + color_len + adrs_len])
        level_id = int(adrs[-data['bit']['q']:], 2)
        block_id = int(adrs[data['bit']['file_id']:-data['bit']['q']], 2)
        color_id = decide_color(line[pre_len: pre_len + color_len], data['channel'][0], data['channel'][1],
                                data['channel'][2])
        for i in range(data['bit']['block'][1]):
            binary_seq[color_id][level_id][block_id] += gc_blck.dna2bit(line[pre_len + color_len + adrs_len + i * data[
                'gc']['block'][0]: pre_len + color_len + adrs_len + (i + 1) * data['gc']['block'][0]])

    # check for empty blocks
    for i in range(3):
        for j in range(len(level_blocks[i])):
            for k in range(level_blocks[i][j]):
                if len(binary_seq[i][j][k]) == 0:
                    print('Some block empty.')

    # concatenation to have list (color) of list (levels)
    for i in range(3):
        for j in range(len(level_blocks[i])):
            for k in range(1, level_blocks[i][j]):
                binary_seq[i][j][0] += binary_seq[i][j][k]
            binary_seq[i][j] = binary_seq[i][j][0]

    print('DNA to levels finished.')
    return binary_seq


def level2bin(huff_dict_inverse_path, binary_seq, file_id):
    """
    Decode binary strings using inverse Huffman dictionary.
    :param huff_dict_inverse_path: str, path of inverse Huffman dictionary
    :param binary_seq: list, encoded binary strings
    :param file_id: int, file index
    :return: binary_seq: list, decoded information
    """
    path = huff_dict_inverse_path + 'file_' + str(file_id) + '_huffman_dictionary_inverse.json'
    dictionary_inverse = json_load(path)
    for i in range(3):
        for j in range(len(binary_seq[i])):
            binary_seq[i][j] = huffman.dec_new(binary_seq[i][j],dictionary_inverse)
    return binary_seq


def dec_img(reads_dna, data):
    """
    Decoding images from DNA consensus reads. Reconstructed images will be automatically stored in folder
    data['path']['output'].
    :param reads_dna: list, consensus reads of DNA oligos
    :param data: dict, general setting read from configuration file
    :return: recon_q: list, reconstructed results
    """
    pre_len = len(data['prefix']['0'])
    suf_len = len(data['suffix']['0'])
    color_len = 3
    adrs_len = data['gc']['id'][0]
    blks_len = data['gc']['block'][0]
    gc_adrs = gc(data['gc']['id'][0], data['gc']['id'][1], data['gc']['id'][2])

    f_list = []
    f_infs = []
    for line in reads_dna:
        adrs = gc_adrs.dna2bit(line[pre_len + color_len: pre_len + color_len + adrs_len])
        file_id = int(adrs[0:data['bit']['file_id']], 2)
        if file_id not in f_list:
            f_list.append(file_id)
            f_infs += [[line]]
        else:
            f_infs[f_list.index(file_id)] += [line]

    recon_q = []
    for file_idx in range(len(f_list)):
        file_id = f_list[file_idx]
        dec_levels = dna2levels(f_infs[file_idx], data)
        dec_bin = level2bin(data['path']['output'], dec_levels, file_id)
        for col_idx in range(3):
            for lv_idx in range(len(dec_bin[col_idx])):
                dec_bin[col_idx][lv_idx] = conversion.diff_dec(dec_bin[col_idx][lv_idx])
        bgr = [[], [], []]
        bgr[0] = conversion.dequantize_new(conversion.level_2_vector_new(dec_bin[0]), data['bit']['q'])
        bgr[1] = conversion.dequantize_new(conversion.level_2_vector_new(dec_bin[1]), data['bit']['q'])
        bgr[2] = conversion.dequantize_new(conversion.level_2_vector_new(dec_bin[2]), data['bit']['q'])
        row = data['size'][str(file_id)][0]
        col = data['size'][str(file_id)][1]
        hilbert_curve = np.load(data['path']['output'] + 'hilbert_curve_' + str(row) + '_' + str(col) + '_dictionary.npy')
        channel_b = np.zeros([row, col])
        channel_g = np.zeros([row, col])
        channel_r = np.zeros([row, col])
        min_itr = min([len(bgr[0]), len(bgr[1]), len(bgr[2]), hilbert_curve.shape[0]])
        for i in range(min_itr):
            channel_b[hilbert_curve[i][0], hilbert_curve[i][1]] = bgr[0][i]
            channel_g[hilbert_curve[i][0], hilbert_curve[i][1]] = bgr[1][i]
            channel_r[hilbert_curve[i][0], hilbert_curve[i][1]] = bgr[2][i]
        recon = [channel_b, channel_g, channel_r]
        recon_q += [recon]
        img_q = cv2.merge(recon)
        cv2.imwrite(data['path']['output'] + str(file_id) + '_recon.png', img_q)
        print('Image reconstruction finished.')

    return recon_q


def run_length(seq):
    """
    Run length analysis of symbols in input string seq.
    """
    max_G_run = 0
    cur_G_run = 0
    max_run = 0
    cur_run = 0
    max_run_symbol = ''
    last_symbol = ''
    for symbol in seq:
        if symbol == 'G':
            cur_G_run += 1
        else:
            if cur_G_run > max_G_run:
                max_G_run = cur_G_run
            cur_G_run = 0
        if symbol == last_symbol:
            cur_run += 1
        else:
            if cur_run > max_run:
                max_run = cur_run
                max_run_symbol = last_symbol
            cur_run = 1
            last_symbol = symbol
    return max_G_run, max_run, max_run_symbol


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Encoding and decoding scheme")
    parser.add_argument("--config", type=str, default="input.json", help="General setting")
    parser.add_argument("--decoding_only", type=bool, default=False, help="Perform only decoding or not")
    parser.add_argument("--read_path", type=str, default="reads_example/seq_without_errors.txt",
                        help="DNA consensus reads if only performing decoding procedure")
    args = parser.parse_args()
    data = json_load(args.config)
    if not args.decoding_only:
        # perform both encoding and decoding procedure
        if not os.path.exists(data['path']['output']):
            # if the output folder does not exist, create the folder
            os.makedirs(data['path']['output'])
        records, reads_dna, ori_q, level_size = enc_dna(data)
        recon_q = dec_img(reads_dna, data)
        print('Encoding-decoding procedure complete.')
        print('=' * 20)

        # sanity check, uncomment the following lines when using conversion.dequantize in dec_img
        # for i in range(len(ori_q)):
        #     if not ((ori_q[i][0] == recon_q[i][0]).all() and (ori_q[i][1] == recon_q[i][1]).all()
        #             and (ori_q[i][2] == recon_q[i][2]).all()):
        #         raise ValueError('Decoded value does not match encoded information.')

        print('Number of oligos in each pool:', level_size)
        # store the correct sequences
        with open(data['path']['output'] + 'seq_without_errors.txt', 'w') as fw:
            for read_line in reads_dna:
                fw.write(read_line + '\n')
        print('Encoded DNA oligos stored.')
    else:
        # if only perform decoding procedure
        if not os.path.exists(data['path']['output'] + args.config):
            # No encoding has been performed yet
            raise ValueError('At least one round of encoding is needed before any decoding procedure.')
        with open(args.read_path, 'r') as fr:
            reads_dna = fr.read().splitlines()
        recon_q = dec_img(reads_dna, data)
        print('Decoding procedure complete.')

