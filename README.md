# cortex_image

## Installation

```
pip install pycortex
pip install nibabel
pip install h5py
git clone https://github.com/Washington-University/HCPpipelines.git
```
[Freesurfer](https://surfer.nmr.mgh.harvard.edu/) must be installed to run mri_convert.

## Usage

Download [NSD](https://naturalscenesdataset.org/) fsaverage betas. 

1. Convert mgh to gii
```
python example_mgh_to_fsaverage.py --nsd_dir /path/to/nsddata_betas/ppdata/ --output_dir /path/to/gii/output/ --subids 1
```

2. Convert gii to fsLR
```
python example_fsaverage_to_fsLR32k.py --fs_dir /path/to/gii/output/ --output_dir /path/to/fsLR/output/ --subids 1
```

3. Convert fsLR to cortex image
```
python example_cortex_to_image.py --output_dir /path/to/cortex/output/ --subids 1 --data_dir /path/to/fsLR/output/ --subids 1
```

