import argparse
import h5py
import os
import re
import subprocess
from pathlib import Path
import nibabel as nib
import numpy as np


# adapted from pycortex
def _make_flatmask(left,right, height=1024):
    from cortex import polyutils
    from PIL import Image, ImageDraw
    pts = np.vstack([left[0], right[0]])
    polys = np.vstack([left[1], right[1]+len(left[0])])

    left, right = polyutils.trace_poly(polyutils.boundary_edges(polys))

    aspect = (height / (pts.max(0) - pts.min(0))[1])
    lpts = (pts[left] - pts.min(0)) * aspect
    rpts = (pts[right] - pts.min(0)) * aspect

    im = Image.new('L', (int(aspect * (pts.max(0) - pts.min(0))[0]), height))
    draw = ImageDraw.Draw(im)
    draw.polygon(lpts[:,:2].ravel().tolist(), fill=255)
    draw.polygon(rpts[:,:2].ravel().tolist(), fill=255)
    extents = np.hstack([pts.min(0), pts.max(0)])[[0,3,1,4]]
    
    return np.array(im).T > 0, extents

def _make_vertex_cache(left,right, height=1024):
    from scipy import sparse
    from scipy.spatial import cKDTree
    flat = np.vstack([left[0], right[0]])
    polys = np.vstack([left[1], right[1]+len(left[0])])
    valid = np.unique(polys)
    fmax, fmin = flat.max(0), flat.min(0)
    size = fmax - fmin
    aspect = size[0] / size[1]
    width = int(aspect * height)
    grid = np.mgrid[fmin[0]:fmax[0]:width*1j, fmin[1]:fmax[1]:height*1j].reshape(2,-1)

    mask, extents = _make_flatmask(left,right, height=height)
    assert mask.shape[0] == width and mask.shape[1] == height

    kdt = cKDTree(flat[valid,:2])
    dist, vert = kdt.query(grid.T[mask.ravel()])
    dataij = (np.ones((len(vert),)), np.array([np.arange(len(vert)), valid[vert]]))
    return sparse.csr_matrix(dataij, shape=(mask.sum(), len(flat)))

def make_flatmap_image(left,right, data, height=1024, recache=False, nanmean=False, **kwargs):

    mask, extents = _make_flatmask(left, right, height=height)
    
    pixmap = _make_vertex_cache(left,right, height=height)

    if data.shape[0] > 1:
        raise ValueError("Input data was not the correct dimensionality - please provide 3D Volume or 2D Vertex data")

    if data.dtype != np.uint8:
        # Convert data to float to avoid image artifacts
        data = data.astype(np.float64)
    if data.dtype == np.uint8:
        img = np.zeros(mask.shape+(4,), dtype=np.uint8)
        img[mask] = pixmap * data.reshape(-1, 4)
        return img.transpose(1,0,2)[::-1], extents
    else:
        badmask = np.array(pixmap.sum(1) > 0).ravel()
        img = (np.nan*np.ones(mask.shape)).astype(data.dtype)
        mimg = (np.nan*np.ones(badmask.shape)).astype(data.dtype)

        # pixmap is a (pixels x voxels) sparse non-negative weight matrix
        # where each row sums to 1

        if not nanmean:
            # pixmap.dot(vec) gives mean of vec across cortical thickness
            mimg[badmask] = pixmap.dot(data.ravel())[badmask].astype(mimg.dtype)
        else:
            # to ignore nans in the weighted mean, nanmean =
            # sum(weights * non-nan values) / sum(weights on non-nan values)
            nonnan_sum = pixmap.dot(np.nan_to_num(data.ravel()))
            weights_on_nonnan = pixmap.dot((~np.isnan(data.ravel())).astype(data.dtype))
            nanmean_data = nonnan_sum / weights_on_nonnan
            mimg[badmask] = nanmean_data[badmask].astype(mimg.dtype)

        img[mask] = mimg

        return img.T[::-1], extents

class GLMtoSurfaces:
    def __init__(self,left,right,height = 1024, resolution='32k',lroi='',rroi='',output_images_dir='surface_images',data_dir='/home1/wangyun/UKB'):
        self.lroi = lroi
        self.rroi = rroi
        self.lroi_data = nib.load(lroi).agg_data()
        self.rroi_data = nib.load(rroi).agg_data()
        self.height = height
        self.output_images_dir = output_images_dir
        self.data_dir = data_dir
        self.left=nib.load(left).agg_data()
        self.right=nib.load(right).agg_data()

        self.left[0][:,0] -= self.left[0].max(0)[0]
        self.right[0][:,0] -= self.right[0].min(0)[0]

        self.lpts, self.lpolys = self.left
        self.rpts, self.rpolys = self.right

    def glm_to_surfaces_pool(self,im_name,temp_left,temp_right):

        if im_name.exists():
            print(f'{im_name} exists')
        else:
            left_data = nib.load(temp_left)
            right_data = nib.load(temp_right)
            len_data=len(left_data.darrays)
            assert len_data==len(right_data.darrays)
            images=[]
            for i in range(len_data):
                temp_left=left_data.darrays[i].data
                temp_right=right_data.darrays[i].data
                if self.lroi!='' and self.rroi!='':
                    temp_left[self.lroi_data==0]=0
                    temp_right[self.rroi_data==0]=0
                lr = np.hstack([temp_left,temp_right]).reshape(1,-1)
                im,ex=make_flatmap_image(self.left,self.right,lr,height=self.height,nanmean=False)
                im[np.isnan(im)] = 0
                if self.lroi!='' and self.rroi!='':
                    idx=np.where(im!=0)
                    im = im[idx[0].min():idx[0].max(),idx[1].min():idx[1].max()]
                im = np.expand_dims(im,axis=0)
                images.append(im)
            
            # normalization
            images=np.concatenate(images,axis=0).astype(np.float32)
            images_mean=np.mean(images,axis=0)
            images_std=np.std(images,axis=0,ddof=1)
            images=(images-images_mean)/images_std
            images[np.isnan(images)]=0
            print(images.shape, images.max(), images.min())

            # print(im_name)
            with h5py.File(f'{im_name}', 'w') as f:
                ds=f.create_dataset('images',images.shape,compression="gzip", compression_opts=9, data=images)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str,required=False,default='/cortex_image/',help='output_dir')
    parser.add_argument('--height', type=int, required=False,default=1024,help='surface image height')
    parser.add_argument('--lroi', type=str, required=False,default='cortex.L.func.gii',help='lroi')
    parser.add_argument('--rroi', type=str, required=False,default='cortex.R.func.gii',help='rroi')
    parser.add_argument('--n_jobs', type=int, required=False, default=1,help='')
    parser.add_argument('--data_dir', type=str,required=False,default='/fsLR32k/',help='fslr data dir')
    parser.add_argument('--left_surf', type=str, required=False,default='S1200.L.flat.32k_fs_LR.surf.gii',help='left_surf')
    parser.add_argument('--right_surf', type=str, required=False,default='S1200.R.flat.32k_fs_LR.surf.gii',help='right_surf')
    parser.add_argument('--subids', type=int, nargs='+',required=False,default=[1],help='subids to convert')

    args = parser.parse_args()
    print(args)

    surf2image = GLMtoSurfaces(height=args.height,lroi=args.lroi,rroi=args.rroi,left=args.left_surf,right=args.right_surf)

    for k,s in enumerate(args.subids):
        print(f'processing subject {s}')
        sub_dir = Path(args.data_dir)/f'subj0{s}'
        sub_out = Path(args.output_dir)/f'subj0{s}'
        sub_out.mkdir(parents=True,exist_ok=True)
        files = []
        for i in range(1,41): # at most 40 sessions
            temp_left=sub_dir/f'l.sess{i:02d}.32k_fs_LR.func.gii'
            temp_right=sub_dir/f'r.sess{i:02d}.32k_fs_LR.func.gii'
            temp_name = sub_out/f'sess{i:02d}.h5'
            if temp_left.exists() and temp_right.exists():
                print(f'processing session {i}')
                surf2image.glm_to_surfaces_pool(temp_name,temp_left, temp_right)