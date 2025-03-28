import os
import nibabel as nib
import numpy as np
from pathlib import Path
import argparse
import re
from tqdm import tqdm
from multiprocessing import Pool

def wrap(fsaverage_name,temp):
    global sub_out,args

    fslr_name = sub_out/f'{temp[0]}.sess{temp[1]}.{args.res}k_fs_LR.func.gii'
    os.system(f'wb_command -metric-resample {fsaverage_name} {Path(args.hcp_dir)}/global/templates/standard_mesh_atlases/resample_fsaverage/fsaverage_std_sphere.{temp[0].upper()}.164k_fsavg_{temp[0].upper()}.surf.gii {Path(args.hcp_dir)}/global/templates/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.{temp[0].upper()}.sphere.{args.res}k_fs_LR.surf.gii ADAP_BARY_AREA {fslr_name} -area-metrics {args.hcp_dir}/global/templates/standard_mesh_atlases/resample_fsaverage/fsaverage.{temp[0].upper()}.midthickness_va_avg.164k_fsavg_{temp[0].upper()}.shape.gii {args.hcp_dir}/global/templates/standard_mesh_atlases/resample_fsaverage/fs_LR.{temp[0].upper()}.midthickness_va_avg.{args.res}k_fs_LR.shape.gii')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--fs_dir', type=str,required=False,default='/fs/',
            help='fsaverage dir')
    parser.add_argument(
            '--hcp_dir', type=str,required=False,default='/HCPpipelines/',
            help='hcp dir')
    parser.add_argument(
            '--output_dir', type=str,required=False,default='/fsLR32k/',
            help='output fslr data dir')
    parser.add_argument(
            '--subids', type=int, nargs='+',required=False,default=[1],
            help='subids to convert')
    parser.add_argument('-r',
            '--res', type=int, required=False, default=32,
            help='resolution, e.g. 32 59 164')
    parser.add_argument(
            '--n_jobs', type=int, required=False, default=64,
            help='')
    args = parser.parse_args()
    print(args)
    
    # load nsd
    for k,s in enumerate(args.subids):
        print(f'processing subject {s}')
        sub_dir = Path(args.fs_dir)/f'subj0{s}'
        sub_out = Path(args.output_dir)/f'subj0{s}'
        sub_out.mkdir(parents=True,exist_ok=True)
        files = list(sub_dir.glob('*gii'))

        for f in tqdm(files):
            temp = re.findall('(?P<lh>[l,r]*)\.sess(?P<ses>[0-9]*)\.func\.gii',f.name)
            wrap(f,temp[0])

