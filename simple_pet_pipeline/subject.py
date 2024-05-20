import numpy as np
import nibabel as nib
import re
import os
import ants
import argparse
import json
from argparse import ArgumentParser
from pathlib import Path
from sys import argv
from glob import glob
from utils import get_file, align, transform, get_tacs
from qc import ImageParam

class Subject():
    def __init__(self, data_dir, out_dir, sub, stx_fn, atlas_fn, labels, clobber=False):
        '''
        Inputs:
            data_dir :  str, directory from which to get pet and mri data for each subject
            out_dir  :  str, directory where outputs will be written
            sub :       str, subject id
            stx_fn :    str, file path to stereotaxic anatomic template
            atlat_fn:   str, file path to stereotaxic atlas
            labels:     dict, labels to use for extracting tacs 
            clobber:    bool, overwrite
        '''
        # Inputs :
        self.data_dir = data_dir
        self.sub = sub
        self.clobber = clobber
        self.roi_labels = labels['roi']
        self.ref_labels = labels['ref']
        self.pet = get_file(data_dir, f'sub-{sub}_*pet.nii.gz')
        self.mri = get_file(data_dir, f'sub-{sub}_*T1w.nii.gz')
        self.stx = stx_fn
        self.atlas_fn = atlas_fn



        # Outputs :
        self.sub_dir = out_dir + os.sep + 'sub-'+ sub
        self.qc_dir = self.sub_dir + os.sep + 'qc/'
        self.pet3d = self.sub_dir + '/' +re.sub('.nii.gz','_3d.nii.gz',os.path.basename(self.pet))
        self.pet_json = get_file(data_dir, f'sub-{sub}_*pet.json')
        self.pet_header = json.load(open(self.pet_json, 'r'))

        self.tacs_csv = self.sub_dir + f'sub-{sub}_TACs.csv'

        self.tacs_qc_plot = self.qc_dir + f'sub-{sub}_TACs.png'

        self.mri2pet_qc_gif = self.qc_dir + f'sub-{sub}_mri2pet.gif'
        self.stx2mri_qc_gif = self.qc_dir + f'sub-{sub}_stx2mri.gif'

        # Class variables 
        self.all_labels = [ i for l in labels.values() for i in l ]
        self.prefix = self.sub_dir + '/' + 'sub-' + sub + '_'



        
        # following variables are defined during <process()>. 
        # not necessary to define these variables here, but helps keeps things clear
        self.atlas_space_pet=None
        self.mri_space_pet=None
        self.pet_space_mri=None
        self.mri2pet_tfm=None
        self.pet2mri_tfm=None
        self.stx_space_mri=None
        self.mri_space_stx=None
        self.stx2mri_tfm=None
        self.mri2stx_tfm=None

        # create output directories
        os.makedirs(self.sub_dir, exist_ok=True)
        os.makedirs(self.qc_dir, exist_ok=True)

        self.frame_duration, self.frame_time_start, self.frame_weight = self.set_frame_times()


    def process(self):
        # Integrate 4D PET image to 3D
        self.pet_to_3d()

        # Align MRI to PET with Rigid Alignment
        self.mri2pet()

        # Align Stereotaxic template to MRI with non-linear SyN transformation
        self.stx2mri()

        # Combine transformations so that we can transform from stereotaxic to PET coord space
        self.stx2pet_tfm = [self.mri2pet_tfm,self.stx2mri_tfm ]

        # Apply stx2pet transformation to stereotaxic atlas
        self.atlas_space_pet = transform(self.prefix, self.pet3d, self.atlas_fn, self.stx2pet_tfm, interpolator='nearestNeighbor', qc_filename=f'{self.qc_dir}/atlas_pet_space.gif', clobber=self.clobber )

        # Apply stx2pet transformation to stereotaxic template
        self.stx_space_pet = transform(self.prefix, self.pet3d, self.stx, self.stx2pet_tfm,  qc_filename=f'{self.qc_dir}/template_pet_space.gif', clobber=self.clobber)

        # Extract time-activity curves (TACs) from PET image using atlas in PET space 
        self.tacs = get_tacs( self.pet, self.atlas_space_pet, self.all_labels ,  self.frame_time_start, self.tacs_csv, self.tacs_qc_plot )

    def set_frame_times(self):
        frame_duration = np.array(self.pet_header['FrameDuration']).astype(float)
        frame_time_start = np.array(self.pet_header['FrameTimesStart']).astype(float)
        frame_weight = frame_duration / np.sum(frame_duration)
        return frame_duration, frame_time_start, frame_weight

    def pet_to_3d(self):
        img = nib.load(self.pet)
        vol = img.get_fdata()
        if len(vol.shape) == 4 : vol = np.sum(vol*self.frame_weight, axis=3)
        nib.Nifti1Image(vol, img.affine).to_filename(self.pet3d)

    def mri2pet(self):
        self.mri_space_pet, self.pet_space_mri, self.mri2pet_tfm, self.pet2mri_tfm = align(self.pet3d, self.mri, transform_method='Rigid', outprefix=f'{self.sub_dir}/sub-{self.sub}_mri2pet_Rigid_', qc_filename = self.mri2pet_qc_gif)
    

    def stx2mri(self):
        self.stx_space_mri, self.mri_space_stx, self.stx2mri_tfm, self.mri2stx_tfm = align(self.mri, self.stx, transform_method='SyN', outprefix=f'{self.sub_dir}/sub-{self.sub}_stx2mri_SyN_', qc_filename = self.stx2mri_qc_gif)
 


