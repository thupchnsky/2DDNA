Programming language: Python 3.6; operating system: Windows 10, CentOS 7.7.1908

The organization of this software is shown below. Detailed usage of each file is included at the begining of the file. 
Note that the raw sequencing files are too large to be included here. The link to those files will be provided if 
requested.

=======first dimension========
Image encoding and decoding: main_encoding_decoding.py
Packages "opencv-python" and "biopython" are needed in this file.
Note that one round of encoding is needed before any kind of decoding, since the Hilbert curve and Huffman dictionary needs
to be generated. The default output folder is "results/". Eight movie posters used in our experiments are included in "input/".
Two sample consensus reads files obtained from sequencing results are included in "reads_example", which can be used for 
direct decoding. The reconstructed images are named as "*_recon.png".

Automatic discoloration detection: auto_detect.py
There are many parameters to tune in this method. Different combinations can be explored to improve the performance. The
combination shown in the code is just an example. Masked images and masks will be stored in different folders. One pair of 
examples is stored in "masked_images/" and "masks/".

Image inpainting: not included here. Please refer to github links:
EdgeConnect: https://github.com/knazeri/edge-connect
GatedConvolution: https://github.com/JiahuiYu/generative_inpainting
Please follow the install instructions provided in the links above to set up these image inpainting methods. Either one of 
them can work with our software, and you can also choose whatever other inpainting methods that you like. When the methods
are correctly configured, put masked images and masks obtained from auto_detect.py into right path. After inpainting, the 
masked regions should be filled in with some reasonable values. One example is stored in "inpainting_results/".

Final smoothing: final_blending.py
Packages "opencv-python" is needed in this file.
There are many parameters to tune in this method. Different combinations can be explored to improve the performance. The
combination shown in the code is just an example. Refined output will be stored in another folder, and adaptive median
filtering can be iteratively applied on the results from last round. One example is stored in "refined_results/".

Image enhancement: not included here. Please refer to github links:
Old Photo Restoration: https://github.com/microsoft/Bringing-Old-Photos-Back-to-Life
Please follow the install instructions provided in the link above to set up this image enhancement algorithm. When the method
is correctly configured, put images obtained from final_blending.py into right path. After enhancement, the blocking effects
coule be completely eliminated. However, this enhancement may not work well on images with granular faces.

=======second dimension========
Nick identification procedure: nick_detect_write.py, nick_detect_erase.py, nick_detect_rewrite.py
Packages "biopython" is needed in these files.
These three files correspond to three rounds of experiments: writing "ILLINOIS", erasing "ILLINOIS" and rewriting "GRAINGER".
There are many parameters to tune in this procedure. Different combinations can be explored to get different results. The
combination shown in the code is just an example, but the conclusion will be the same as the one shown in our paper.