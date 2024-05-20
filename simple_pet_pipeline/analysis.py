import numpy as np
import nibabel as nib
import pandas as pd
import re
import os
import ants
import argparse
from argparse import ArgumentParser
from pathlib import Path
from sys import argv
from glob import glob

def get_stats_for_labels(vol, atlas, labels):
    total=0
    n=0
    maximum=np.min(vol)

    for l in labels :
        idx = atlas == l
        label_n = np.sum(idx)
        if label_n > 0 :
            pet_values_in_label = vol[ idx ]

            total += np.sum(vol[idx])
            max_in_label = np.max(vol[idx])
            maximum = max_in_label if max_in_label > maximum else maximum

            n += label_n
        else :
            print(f'Error: label {l} not found in atlas volume where it was expected. Skipping')
            exit(1)

    average = total / n
    return average, maximum

def tumor_striatum_analysis(subject, roi_labels, ref_labels):
    print('\tTumor Striatum Analysis')
    print(f'\t\tPET:\t{subject.pet3d}')
    print(f'\t\tAtlas:\t{subject.atlas_space_pet}')

    
    subject.pet_suvr = subject.prefix+'suvr_max.nii.gz'
    subject.suvr_csv = subject.prefix+'suvr_max.json'

    if not os.path.exists(subject.pet_suvr) or not os.path.exists(subject.suvr_csv) or subject.clobber :
        pet_hd = nib.load(subject.pet3d)
        pet_vol = pet_hd.get_fdata()


        atlas_hd = nib.load(subject.atlas_space_pet)
        atlas_vol = np.rint(atlas_hd.get_fdata()).astype(int)

        roi_avg, roi_max = get_stats_for_labels(pet_vol, atlas_vol, roi_labels)
        ref_avg, ref_max = get_stats_for_labels(pet_vol, atlas_vol, ref_labels)

        suvr_max = pet_vol / ref_max

        tumor_labels = suvr_max > 1

        tumor_avg, tumor_max = get_stats_for_labels(pet_vol, atlas_vol, tumor_labels)
        
        ts_ratio = np.round(tumor_max / roi_max,3)

        nib.Nifti1Image(suvr_max, nib.load(subject.pet3d).affine).to_filename(subject.pet_suvr)
        
        suvr_dict = {'sub':[subject.sub],'tumor_max':[tumor_max],'tumor_avg':[tumor_avg], 'striatum_max':[roi_max], 'straitum_avg':[roi_avg], 'ts_ratio':[ts_ratio]}
        subject.suvr_df = pd.DataFrame(suvr_dict)
        subject.suvr_df.to_csv(subject.suvr_csv, index=False)

        #print(f'\t\tTumor   \tAvg:\t{tumor_avg:.5}\tMax:\t{tumor_max:.5}')
        #print(f'\t\tStriatum\tAvg:\t{roi_avg:.5}\tMax:\t{roi_max:.5}')
        #print(f'\t\tTumor_max / Striatum_max : {ts_ratio:.2}')
        print(f'\t\tWriting SUVR_max volume to {subject.pet_suvr}')
    else :
        subject.suvr_df = pd.read_csv(subject.suvr_csv)

    return subject




