#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This file is for final refinement after image inpainting, including bilateral filtering and adaptive median
        filtering.
        fig_path: str, path of target image;
        mask_path: str, path of the mask image;
        save_flag: bool, whether to save the refined image or not;
        adaptive: bool, whether to perform median filtering over previous refinement or on inapinted result without
            refinement. Please use images in folder refined_images/ if adaptive=True, since bilateral filtering will
            not be performed in this case;
        fig_out_path: str, output path for refined images

    Usage:
        if dealing with raw inpainted results, use:
            python final_blending.py --fig_path="inpainting_results/4_inpaint.png" --mask_path="masks/mask_4_recon.png"
        if dealing with previous round of refinement, then use:
            python final_blending.py --fig_path="inpainting_results/4_inpaint.png" --mask_path="masks/mask_4_recon.png"
"""

import cv2
from matplotlib import pyplot as plt
from PIL import Image
import numpy as np
import copy
import os
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Blending procedure")
    parser.add_argument("--fig_path", type=str, default="inpainting_results/4_inpaint.png",
                        help="Path of target image")
    parser.add_argument("--mask_path", type=str, default="masks/mask_4_recon.png", help="Path of mask image")
    parser.add_argument("--save_flag", type=bool, default=True, help="Save the results or not")
    parser.add_argument("--adaptive", type=bool, default=False,
                        help="Perform median filtering over previous refinement or not")
    parser.add_argument("--fig_out_path", type=str, default="refined_results/", help="Output path for refined images")
    args = parser.parse_args()
    # Load mask image
    mask = np.array(Image.open(args.mask_path))
    m, n = mask.shape
    # Load inpainted image
    img_edgecon = np.array(Image.open(args.fig_path))
    if not args.adaptive:
        # If work on inpainting results first
        # Bilateral filtering. This set of coefficients (9,30,30) works well for 4_inpaint.png
        img = cv2.bilateralFilter(img_edgecon, 9, 30, 30)
    else:
        # If work on previous round of refinement, skip bilateral filtering this time
        img = copy.copy(img_edgecon)
    # Median filtering around masks
    for i in range(m):
        for j in range(n):
            if mask[i, j] == 255:
                i_start = min(i + 1, m - 1)
                j_start = min(j + 1, n - 1)
                lst = []
                while i_start < m and j_start < n:
                    lst = list(mask[i_start, j:j_start + 1])
                    if sum(lst) > 0:
                        i_start += 1
                        continue
                    else:
                        lst = list(mask[i:i_start+1, j_start])
                        if sum(lst) > 0:
                            j_start += 1
                            continue
                    break
                tmp = mask[i-1:i_start+1, j-1:j_start+1]
                # the size of the found mask can be adjusted for different image to obtain better performance
                if tmp.size > 25:
                    img[i:i_start, j:j_start, :] = cv2.medianBlur(img[i:i_start, j:j_start, :], 3)
                    mask[i:i_start, j:j_start] = 0
                # if tmp.size > 50:
                #     img[i-1:i_start+1, j-1:j_start+1, :] = cv2.medianBlur(img[i-1:i_start+1, j-1:j_start+1, :], 3)
                #     mask[i-1:i_start+1, j-1:j_start+1] = 0

    # Uncomment the following lines to see what the refinement looks like, comparing with input image
    plt.subplot(121), plt.imshow(img_edgecon), plt.title('Before Processing')
    plt.xticks([]), plt.yticks([])
    plt.subplot(122), plt.imshow(img), plt.title('After Processing')
    plt.xticks([]), plt.yticks([])
    plt.show()

    if args.save_flag:
        fig_name = args.fig_path.split('/')[-1]
        if not os.path.exists(args.fig_out_path):
            # If the output folder does not exist, create the folder
            os.makedirs(args.fig_out_path)
        o1 = Image.fromarray(img.astype(np.uint8))
        if not args.adaptive:
            o1.save(args.fig_out_path + 'refined_' + fig_name)
        else:
            o1.save(args.fig_path)
