import os
import nibabel as nib
import numpy as np
from pathlib import Path
import argparse
import re


def wrap_tobiggii(f):
    global sub_out,args

    temp = re.findall('(?P<l>[l,r])h.betas_session(?P<s>[0-9]*).mgh',f.name)

    if len(temp)>0:
        temp = temp[0]
        gii_file = sub_out/f'{temp[0]}.sess{temp[1]}.func.gii'
        os.system(f'mri_convert {f} {gii_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--nsd_dir', type=str,required=False,default='/natural-scenes-dataset/nsddata_betas/ppdata/',
            help='nsd betas ppdata dir')
    parser.add_argument(
            '--hcp_dir', type=str,required=False,default='/HCPpipelines/',
            help='hcp pipeline dir')
    parser.add_argument(
            '--output_dir', type=str,required=False,default='/fs/',
            help='output fsaverage gii data dir')
    parser.add_argument(
            '--subids', type=int, nargs='+',required=False,default=[1],
            help='subids to convert')
    parser.add_argument('-r',
            '--res', type=int, required=False, default=32,
            help='resolution, e.g. 32 59 164')
    parser.add_argument(
            '--n_jobs', type=int, required=False, default=16,
            help='')
    args = parser.parse_args()
    print(args)

    for k,s in enumerate(args.subids):
        print(f'processing subject {s}')
        sub_dir = Path(args.nsd_dir)/f'subj0{s}'
        fs_dir = sub_dir/'fsaverage/betas_fithrf_GLMdenoise_RR/'
        sub_out = Path(args.output_dir)/f'subj0{s}'
        sub_out.mkdir(parents=True,exist_ok=True)
        files = list(fs_dir.glob('*betas*'))
        for f in files:
            wrap_tobiggii(f)
