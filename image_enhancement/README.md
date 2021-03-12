# Old Photo Restoration

Please refer to the original github link for more information: https://github.com/microsoft/Bringing-Old-Photos-Back-to-Life.

## Requirement
The code is tested on Ubuntu with Nvidia GPUs and CUDA installed. Python>=3.6 is required to run the code.

## Installation

Clone the Synchronized-BatchNorm-PyTorch repository for

```
cd Face_Enhancement/models/networks/
git clone https://github.com/vacancy/Synchronized-BatchNorm-PyTorch
cp -rf Synchronized-BatchNorm-PyTorch/sync_batchnorm .
cd ../../../
```

```
cd Global/detection_models
git clone https://github.com/vacancy/Synchronized-BatchNorm-PyTorch
cp -rf Synchronized-BatchNorm-PyTorch/sync_batchnorm .
cd ../../
```

Download the landmark detection pretrained model

```
cd Face_Detection/
wget http://dlib.net/files/shape_predictor_68_face_landmarks.dat.bz2
bzip2 -d shape_predictor_68_face_landmarks.dat.bz2
cd ../
```

Download the pretrained model from Azure Blob Storage, put the file `Face_Enhancement/checkpoints.zip` under `./Face_Enhancement`, and put the file `Global/checkpoints.zip` under `./Global`. Then unzip them respectively.

```
cd Face_Enhancement/
wget https://facevc.blob.core.windows.net/zhanbo/old_photo/pretrain/Face_Enhancement/checkpoints.zip
unzip checkpoints.zip
cd ../
cd Global/
wget https://facevc.blob.core.windows.net/zhanbo/old_photo/pretrain/Global/checkpoints.zip
unzip checkpoints.zip
cd ../
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## How to use?

**Note**: GPU can be set 0 or 0,1,2 or 0,2; use -1 for CPU

### 1) Full Pipeline

You could easily restore the old photos with one simple command after installation and downloading the pretrained model.

For images without scratches:

```
python run.py --input_folder [test_image_folder_path] \
              --output_folder [output_path] \
              --GPU 0
```

For scratched images:

```
python run.py --input_folder [test_image_folder_path] \
              --output_folder [output_path] \
              --GPU 0 \
              --with_scratch
```

Note: Please try to use the absolute path. The final results will be saved in `./output_path/final_output/`. You could also check the produced results of different steps in `output_path`.

### 2) Scratch Detection

Currently we don't plan to release the scratched old photos dataset with labels directly. If you want to get the paired data, you could use our pretrained model to test the collected images to obtain the labels.

```
cd Global/
python detection.py --test_path [test_image_folder_path] \
                    --output_dir [output_path] \
                    --input_size [resize_256|full_size|scale_256]
```

### 3) Global Restoration

A triplet domain translation network is proposed to solve both structured degradation and unstructured degradation of old photos.

```
cd Global/
python test.py --Scratch_and_Quality_restore \
               --test_input [test_image_folder_path] \
               --test_mask [corresponding mask] \
               --outputs_dir [output_path]

python test.py --Quality_restore \
 --test_input [test_image_folder_path] \
 --outputs_dir [output_path]
```

### 4) Face Enhancement

We use a progressive generator to refine the face regions of old photos. More details could be found in our journal submission and `./Face_Enhancement` folder.

## Citation

If you use this component, please consider citing the following papers

```
@inproceedings{wan2020bringing,
title={Bringing Old Photos Back to Life},
author={Wan, Ziyu and Zhang, Bo and Chen, Dongdong and Zhang, Pan and Chen, Dong and Liao, Jing and Wen, Fang},
booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition},
pages={2747--2757},
year={2020}
}
```

```
@misc{2009.07047,
Author = {Ziyu Wan and Bo Zhang and Dongdong Chen and Pan Zhang and Dong Chen and Jing Liao and Fang Wen},
Title = {Old Photo Restoration via Deep Latent Space Translation},
Year = {2020},
Eprint = {arXiv:2009.07047},
}
```