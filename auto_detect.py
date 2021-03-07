#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for automatic discoloration detection. fig_path: str, path of target image;
        save_flag: bool, whether to save the image and mask or not;
        mask_expand_flag: bool, whether to perform expansion on initial mask or not;
        fig_out_path: str, output path for masked images;
        mask_out_path: str, Output path for masks

    Usage:
        python auto_detect.py --fig_path="results/4_recon.png"
"""

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import copy
import os
import argparse

# num_detect_01, num_detect_02 and num_detect_12 define the number of bins that we assume caused by errors in the
# histogram. Change these three values can affect performance on different images. The larger the value is, the more
# pixels will be masked.
# block_offset and pre_thres is to expand the mask after initial computation to cover other potential defects
# The following setting works well for 4_recon.png
num_detect_01 = 15
num_detect_02 = 15
num_detect_12 = 15
block_offset = 5
per_thres = 0.6


def hist_detect(a, b, num_detect):
    diff = a - b
    n1, inter1, _ = plt.hist(diff.flatten(), bins=30)
    inter_bound = np.zeros((num_detect, 2))
    o1 = np.argsort(n1)
    for i in range(num_detect):
        inter_bound[i, 0] = inter1[o1[i]]
        inter_bound[i, 1] = inter1[o1[i] + 1]
    row_idx = []
    col_idx = []
    for i in range(num_detect):
        tmp = np.where((diff >= inter_bound[i, 0]) & (diff < inter_bound[i, 1]))
        row_idx += list(tmp[0])
        col_idx += list(tmp[1])
    return row_idx, col_idx


def err_pos_combine_naive(r1, r2, r3, c1, c2, c3):
    # naive combination of erroneous pixel positions
    return r1 + r2 + r3, c1 + c2 + c3


def err_pos_combine_majorvote(r1, r2, r3, c1, c2, c3):
    # combination of erroneous pixel positions using majority vote strategy
    r_g = set([(i, j) for i, j in zip(r1, c1)])
    r_b = set([(i, j) for i, j in zip(r2, c2)])
    g_b = set([(i, j) for i, j in zip(r3, c3)])
    # r_err = (r_g.intersection(r_b)).difference(g_b)
    # g_err = (r_g.intersection(g_b)).difference(r_b)
    # b_err = (r_b.intersection(g_b)).difference(r_g)
    r_err = r_g.intersection(r_b)
    g_err = r_g.intersection(g_b)
    b_err = r_b.intersection(g_b)
    err_set = r_err.union(g_err, b_err)
    r = []
    c = []
    for err_pos in err_set:
        r.append(err_pos[0])
        c.append(err_pos[1])
    return r, c


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automatic discoloration detection")
    parser.add_argument("--fig_path", type=str, default="results/4_recon.png", help="Path of target image")
    parser.add_argument("--save_flag", type=bool, default=True, help="Save the results or not")
    parser.add_argument("--mask_expand_flag", type=bool, default=True, help="Perform expansion on initial mask or not")
    parser.add_argument("--fig_out_path", type=str, default="masked_images/", help="Output path for masked images")
    parser.add_argument("--mask_out_path", type=str, default="masks/", help="Output path for masks")
    args = parser.parse_args()
    # Load the image
    img_raw = np.array(Image.open(args.fig_path))

    # Uncomment the following lines to see what each channel looks like
    # plt.subplot(221), plt.imshow(img_raw[:, :, 0]), plt.title('R Channel')
    # plt.xticks([]), plt.yticks([])
    # plt.subplot(222), plt.imshow(img_raw[:, :, 1]), plt.title('G Channel')
    # plt.xticks([]), plt.yticks([])
    # plt.subplot(223), plt.imshow(img_raw[:, :, 2]), plt.title('B Channel')
    # plt.xticks([]), plt.yticks([])
    # plt.subplot(224), plt.imshow(img_raw), plt.title('Original Image')
    # plt.xticks([]), plt.yticks([])
    # plt.show()

    img = img_raw
    # Find pixel positions with lowest value frequencies
    r1, c1 = hist_detect(img[:, :, 0], img[:, :, 1], num_detect_01)
    r2, c2 = hist_detect(img[:, :, 0], img[:, :, 2], num_detect_02)
    r3, c3 = hist_detect(img[:, :, 1], img[:, :, 2], num_detect_12)
    # Use either naive or majorvote strategy
    r, c = err_pos_combine_naive(r1, r2, r3, c1, c2, c3)
    # r, c = err_pos_combine_majorvote(r1, r2, r3, c1, c2, c3)

    # Initialize the mask
    mask = np.zeros((img.shape[0], img.shape[1]))
    img_copy = copy.copy(img_raw)
    for i in range(len(r)):
        mask[r[i], c[i]] = 1

    mask_new = copy.copy(mask)
    if args.mask_expand_flag:
        for i in range(0, img.shape[0], block_offset):
            for j in range(0, img.shape[1], block_offset):
                if np.sum(mask[i:min(img.shape[0], i+block_offset),
                               j:min(img.shape[1], j+block_offset)]) >= (block_offset**2 * per_thres):
                    mask_new[max(0, round(i-block_offset/2)): min(img.shape[0], round(i+block_offset)),
                             max(0, round(j-block_offset/2)):min(img.shape[1], round(j+block_offset))] = 1

    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if mask_new[i, j] == 1:
                img_copy[i, j, :] = 255
                mask_new[i, j] = 255
            elif mask_new[i, j] != 0:
                print('Value error')

    # Uncomment the following lines to see what masked image looks like, comparing with original image
    # plt.subplot(121), plt.imshow(img), plt.title('With Error')
    # plt.xticks([]), plt.yticks([])
    # plt.subplot(122), plt.imshow(img_copy), plt.title('Error Positions Masked as White')
    # plt.xticks([]), plt.yticks([])
    # plt.show()

    if args.save_flag:
        fig_name = args.fig_path.split('/')[-1]
        if not os.path.exists(args.fig_out_path):
            # If the output folder does not exist, create the folder
            os.makedirs(args.fig_out_path)
        if not os.path.exists(args.mask_out_path):
            # If the output folder does not exist, create the folder
            os.makedirs(args.mask_out_path)
        # Save masked image
        o1 = Image.fromarray(img_copy.astype(np.uint8))
        o1.save(args.fig_out_path + 'masked_' + fig_name)
        # Save mask
        o2 = Image.fromarray(mask_new.astype(np.uint8), 'L')
        o2.save(args.mask_out_path + 'mask_' + fig_name)

