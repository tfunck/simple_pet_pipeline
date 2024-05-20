import numpy as np
import nibabel as nib
import re
import os
import ants
import argparse
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from sys import argv
from glob import glob
from subject import Subject
from analysis import tumor_striatum_analysis
from utils import get_file



def find_subject_ids(data_dir):
    get_id = lambda fn: re.sub('sub-','',os.path.basename(fn).split('_')[0])
    pet_images_list = Path(data_dir).rglob('*_pet.nii.gz')
    return [ get_id(fn) for fn in pet_images_list ]

        
def get_parser():
    parser = ArgumentParser(usage="useage: ")
    parser.add_argument("-i",dest="data_dir", default='pediatric/', help="Path for input file directory")
    parser.add_argument("-o",dest="out_dir", default='output/', help="Path for output file directory")
    parser.add_argument("-s",dest="stx_fn", default='atlas/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz', help="Path for stereotaxic template file")
    parser.add_argument("-a",dest="atlas_fn", default='atlas/dka.nii.gz', help="Path for stereotaxic label file")
    parser.add_argument("--ref-labels",dest="ref_labels", type=int, help="Label values to use for Quantification", default=[2], nargs='+' )
    parser.add_argument("--roi-labels",dest="roi_labels", type=int, help="Label values to use for Quantification", default=[11,12,13], nargs='+' )
    return parser


if __name__ == '__main__' :

    
    opts = get_parser().parse_args()
    print('\n Pediatric FDOPA Pipeline\n ________________________\n')
    print('\tOptions')
    print('\t\tData directory:', opts.data_dir)
    print('\t\tOutput directory:', opts.out_dir)
    print('\t\tTemplate:',opts.stx_fn)
    print('\t\tAtlas:', opts.atlas_fn)
    print(f'\t\tRef. Labels: {opts.ref_labels}')
    print(f'\t\tStr. Labels: {opts.roi_labels}')
    print()

    tumor_striatum_csv = opts.out_dir+os.sep+'tumor_striatum.csv'
    print('\tOutputs')
    print('\t\tTumor striatum ratio csv:'+tumor_striatum_csv)
    print()

    subject_id_list = find_subject_ids(opts.data_dir)
    print('\tRuntime parameters:')
    print(f'\t\tSubject IDs: {subject_id_list}')
    print()

    labels = {'roi':opts.roi_labels, 'ref': opts.ref_labels}

    # Create a list of instances of the Subject class. 
    subject_list = [ Subject(opts.data_dir, opts.out_dir, sub, opts.stx_fn, opts.atlas_fn, labels) for sub in subject_id_list ]
    
    # Do initial processing for each subject (e.g., alignment to MRI and stereotaxic atlas)
    [ subj.process() for subj in subject_list ] 

    # Do analysis to find maximum tumor and striatum PET values
    subject_list = [ tumor_striatum_analysis(subj, opts.roi_labels, opts.ref_labels) for subj in subject_list ]
    tumor_striatum_df = pd.concat([ subject.suvr_df for subject in subject_list ])
    tumor_striatum_df.to_csv(tumor_striatum_csv, index=False)

    print()

