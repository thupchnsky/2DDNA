# Official implementation of 2DDNA image storage platform

### [Paper (Under Review)](https://www.biorxiv.org/content/10.1101/2021.02.22.432304v2)

Programming language: Python 3.6. Tested on operating systems: Windows 10, CentOS 7.7.1908

## Organization
The organization of this software is shown below. Detailed usage of each file is included at the begining of the file. 
Note that the raw sequencing files are too large to be included here. The link to those files will be provided if 
requested.

## First Dimension (Image Content)
Image encoding and decoding: `main_encoding_decoding.py`. Packages `opencv-python` and `biopython` are needed. Note that at least one round of encoding is needed before any kind of decoding, since the Hilbert curve and Huffman dictionary needs to be generated. The default output folder is `results/`. Eight movie posters used in our experiments are included in `input/`. Two sample consensus reads files obtained from sequencing results are included in `reads_example/` (one with errors and one without), which can be used for direct decoding. The reconstructed images are named as `*_recon.png`.

## Automatic Discoloration Detection
Automatic discoloration detection: `auto_detect.py`. There are many parameters to tune in this method. Different combinations can be explored to improve the performance. The
combination shown in the code is just an example. Masked images and masks will be stored in different folders. One pair of examples is stored in `masked_images/` and `masks/`.

## Image Inpainting
There are many off-the-shelf and well developed image inapinting methods avaliable online. Two methods that we adopted into our framework can be found here:<br>
<strong>EdgeConnect</strong>: https://github.com/knazeri/edge-connect<br>
<strong>GatedConvolution</strong>: https://github.com/JiahuiYu/generative_inpainting<br>
Please follow the install instructions provided in the links above to set up these image inpainting methods. Either one of them can work well with our software, and you can also choose whatever other inpainting methods that you like. When the methods are correctly configured, put masked images and masks obtained from `auto_detect.py` into right path. After inpainting, the masked regions should be filled with some reasonable values. One example is stored in `inpainting_results/`.

## Final Smoothing
Final smoothing: `final_blending.py`. Package `opencv-python` is needed in this file. There are many parameters to tune in this method. Different combinations can be explored to improve the performance. The combination shown in the code is just an example. Refined output will be stored in another folder, and adaptive median filtering can be iteratively applied on the results from last round. One example is stored in `refined_results/`.

## Image Enhancement
For image enhancement component please refer to following github link:<br>
<strong>Old Photo Restoration</strong>: https://github.com/microsoft/Bringing-Old-Photos-Back-to-Life<br>
Please follow the install instructions provided in the link above to set up this image enhancement algorithm. When the method is correctly configured, put images obtained from `final_blending.py` into right path. After enhancement, the blocking effects coule be completely eliminated. However, this enhancement may <strong>not</strong> work well on images with granular faces.

## Second Dimension (Image Metadata)
Nick identification procedure: `nick_detect_write.py`, `nick_detect_erase.py`, `nick_detect_rewrite.py`. Package "biopython" is needed in these files. These three files correspond to three rounds of experiments: writing "ILLINOIS", erasing "ILLINOIS" and rewriting "GRAINGER". There are many parameters to tune in this procedure. Different combinations can be explored to get different results. The combination shown in the code is just an example, but the conclusion will be the same as the one shown in our paper.

## Contact
Please contact Chao Pan (chaopan2@illinois.edu) if you have any questions.
