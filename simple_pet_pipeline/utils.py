import matplotlib 
matplotlib.rcParams['figure.facecolor'] = '1.'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import re
import os
import ants
import argparse
import pandas as pd
import seaborn as sns
from qc import ImageParam
from argparse import ArgumentParser
from pathlib import Path
from sys import argv
from glob import glob

### Utility functions

def get_tacs(pet, atlas, target_labels, times, tac_csv, qc_png=None, clobber=False):
    df = pd.DataFrame({'frame':[], 'label':[],'region':[], 'value':[]})
    pet_vol = nib.load(pet).get_fdata()
    atlas_vol = np.rint(nib.load(atlas).get_fdata()).astype(int)


    print('\tExtract Time Activity Curves')
    print('\t\tPet:',pet)
    print('\t\tAtlas:',atlas)
    print('\t\tLabels:', target_labels)
    print('\t\tTAC csv:', tac_csv)
    print('\t\tTAC QC:', qc_png)
    print()
    clobber=True
    if not os.path.exists(tac_csv) or clobber :

        rows = []
        
        for label in target_labels:
            for frame in range(pet_vol.shape[3]) :

                pet_frame = pet_vol[:,:,:,frame]

                assert np.sum(atlas_vol==label) > 0, f'Error: could not find {label} in atlas'

                values = pet_frame[atlas_vol==label]

                if np.sum(np.abs(values)) > 0 :
                    avg = np.mean(values)
                else:
                    avg=0

                new_row = pd.DataFrame({'time':[times[frame]],'frame':[frame], 'region':[label], 'value':[avg]})
                rows.append(new_row)

        df = pd.concat(rows)

        if type(qc_png) == str :
            plt.title('Time Activity Curves')
            sns.lineplot(data=df, x='time', y='value', hue='region', marker=True, dashes=False, palette=sns.color_palette("husl", len(target_labels))) 
            plt.savefig(qc_png)

        df.to_csv(tac_csv)

    else :
        df = pd.read_csv(tac_csv)

    return df

def get_file(data_dir, string):
    
    lst = list(Path(data_dir).rglob(string))
    if len(lst)==1 :
        return str(lst[0])
    else :
        print('Could not find single file for string')
        print(lst)
        exit(1)


def transform(prefix, fx, mv, tfm, interpolator='linear', qc_filename=None, clobber=False):
    print('\tTransforming')
    print('\t\tFixed',fx)
    print('\t\tMoving',mv)
    print('\t\tTransformations:', tfm)
    print('\t\tQC:', qc_filename)
    print()
    out_fn = prefix +  re.sub('.nii.gz','_rsl.nii.gz', os.path.basename(mv))
    if not os.path.exists(out_fn) or clobber :
        img_rsl = ants.apply_transforms(fixed=ants.image_read(fx), 
                                        moving=ants.image_read(mv), 
                                        transformlist=tfm,
                                        interpolator=interpolator,
                                        verbose=True
                                        )
        ants.image_write( img_rsl, out_fn )
        
        if type(qc_filename) == str :
            ImageParam(fx, qc_filename, out_fn, duration=600,  nframes=15, dpi=400, alpha=[0.4]   ).volume2gif()
    return out_fn


def align(fx, mv, transform_method='SyN', init=[], outprefix='', qc_filename=None) :
   
    warpedmovout =  outprefix + 'fwd.nii.gz'
    warpedfixout =  outprefix + 'inv.nii.gz'
    fwdtransforms = outprefix+'Composite.h5'
    invtransforms = outprefix+'InverseComposite.h5'

    print(f'\tAligning\n\t\tFixed: {fx}\n\t\tMoving: {mv}\n\t\tTransform: {transform_method}')
    print(f'\t\tQC: {qc_filename}\n')
    output_files = warpedmovout, warpedfixout, fwdtransforms, invtransforms
    if False in [os.path.exists(fn) for fn in output_files ] :
        out = ants.registration(fixed = ants.image_read(fx), 
                                moving = ants.image_read(mv), 
                                type_of_transform = transform_method, 
                                init=init,
                                verbose=True,
                                outprefix=outprefix,
                                write_composite_transform=True
                                )
        ants.image_write(out['warpedmovout'], warpedmovout)
        ants.image_write(out['warpedfixout'], warpedfixout)
        
        if type(qc_filename) == str :
            ImageParam(fx, qc_filename, warpedmovout, duration=600,  nframes=15, dpi=400, alpha=[0.3],  edge_2=1, cmap1=plt.cm.Greys, cmap2=plt.cm.Reds ).volume2gif()

    return output_files
