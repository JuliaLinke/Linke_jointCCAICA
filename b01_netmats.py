#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 10:36:20 2019

@author: winkleram
"""
import nibabel
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import imageio

rootdir = os.path.join('/','data','EDB','MErest')
dataset = "NIMH"
partial = True
remove_vols = 0
fd_thr = .5

def loadGifti(fname, NonSteadyState=0, icres=7): # ============================
    gii = nibabel.load(fname)
    gii_data = [d.data[:,None] for d in gii.darrays]
    gii_data = np.concatenate(gii_data, axis=1).T
    nV = 4**icres*10 + 2
    gii_data = gii_data[:,0:nV]
    return gii_data[NonSteadyState:,:]

def loadNifti(fname, NonSteadyState=0): # =====================================
    n = nibabel.load(fname)
    naff   = n.affine
    img4d  = n.get_fdata();
    imgsiz = img4d.shape
    if len(imgsiz) == 4:
        img4d  = img4d[:,:,:,NonSteadyState:]
        imgsiz = img4d.shape
        img2d  = np.reshape(img4d, (np.prod(imgsiz[0:3]), imgsiz[-1]), order='F').T
    else:
        img2d  = np.reshape(img4d, (np.prod(imgsiz[0:3]), 1), order='F').T
    return img2d, imgsiz, naff

def unWrap(netmat, side='lower'): # ===========================================
    if side.lower() == 'lower':
        uw  = netmat[np.tril_indices(netmat.shape[0], k=-1, m=netmat.shape[1])]
    elif side.lower() == 'upper':
        uw  = netmat[np.triu_indices(netmat.shape[0], k= 1, m=netmat.shape[1])]
    elif side.lower() == 'both':
        uwl = netmat[np.tril_indices(netmat.shape[0], k=-1, m=netmat.shape[1])]
        uwu = netmat[np.triu_indices(netmat.shape[0], k= 1, m=netmat.shape[1])]
        uw  = np.concatenate((uwl, uwu))
    return uw

# ======== [ Main ] ===========================================================
if __name__ == '__main__':
    
    listSubjects = [d for d in os.listdir(os.path.join(rootdir,'derivatives_{}'.format(dataset),'fmriprep')) if (os.path.isdir(os.path.join(rootdir,'derivatives_{}'.format(dataset),'fmriprep',d)) and d[0:4] == 'sub-')]
    listSubjects.sort()
    fmriprepdir  = os.path.join(rootdir, 'derivatives_{}'.format(dataset), 'fmriprep')
    netmatsdir   = os.path.join(fmriprepdir, '..', 'netmats_tmp')
    
    # List of subcortical structures of interest
    # 26/58 = Accumbens
    # 18/54 = Amygdala
    # 11/50 = Caudate
    # 17/53 = Hippocampus
    # 13/52 = Pallidum
    # 12/51 = Putamen
    # 10/49 = Thalamus
    # 28/60 = Ventral DC
    aseg_list = {
            'L': [26, 18, 11, 17, 13, 12, 10, 28],
            'R': [58, 54, 50, 53, 52, 51, 49, 60]}
    
    # Load the parcellation "annot" files in fsaverage space
    annot = {}
    ctab  = {}
    names = {}
    [annot['L'], ctab['L'], names['L']] = nibabel.freesurfer.io.read_annot(os.path.join(rootdir ,'atlases', 'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label','lh.Schaefer2018_200Parcels_17Networks_order.annot'))
    [annot['R'], ctab['R'], names['R']] = nibabel.freesurfer.io.read_annot(os.path.join(rootdir ,'atlases', 'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label','rh.Schaefer2018_200Parcels_17Networks_order.annot'))
    
    # For each subject
    for sidx, subj in enumerate(listSubjects[175:]):
        print(sidx, subj)
        
        # For each session (if existing)
        sesList = glob.glob(os.path.join(fmriprepdir, subj, 'ses-*'))
        if not sesList:
            sesList = ['']
        else:
            for idx, ses in enumerate(sesList):
                sesList[idx] = os.path.split(sesList[idx])[1]
        for ses in sesList:
            
            # For each run (if existing)
            runList = glob.glob(os.path.join(fmriprepdir, subj, ses, 'func', '{}{}_task-rest_run-*_desc-confounds_regressors.tsv'.format(subj, ses if not ses else '_{}'.format(ses))))
            if not runList:
                runList = ['']
            else:
                for idx, run in enumerate(runList):
                    tmp = os.path.split(runList[idx])[1]
                    runList[idx] = tmp.split('_')[2]
            for run in runList:
                # Load the confounds
                confounds = pd.read_csv(os.path.join(fmriprepdir, subj, ses, 'func', '{}{}_task-rest{}_desc-confounds_regressors.tsv'.format(subj, ses if not ses else '_{}'.format(ses), run if not run else '_{}'.format(run))), delimiter='\t')
                confounds['std_dvars'][0]              = 0
                confounds['dvars'][0]                  = 0
                confounds['framewise_displacement'][0] = 0
                Z  = np.concatenate((np.array(confounds.loc[:, confounds.columns != 'global_signal']), np.ones((confounds.shape[0],1))), axis=1)
                
                if remove_vols > 0:
                    Z = Z[remove_vols:,:]
                    confounds.iloc[remove_vols:]
                
                if fd_thr > 0:
                    idx_tp = np.array(confounds['framewise_displacement'] < fd_thr)
                    n_tp   = sum(idx_tp)
                    Z = Z[idx_tp,:]

                # For each hemisphere
                surf_func = {}
                surf_parc = {}
                for hemi in ['L', 'R']:
                    
                    # Load functional data in fsaverage space (surface)
                    surf_func[hemi] = loadGifti(os.path.join(fmriprepdir, subj, ses, 'func', '{}{}_task-rest{}_space-fsaverage_hemi-{}.func.gii'.format(subj, ses if not ses else '_{}'.format(ses), run if not run else '_{}'.format(run), hemi)))
        
                    if remove_vols > 0:
                        surf_func[hemi] = surf_func[hemi][remove_vols:,:]
        
                    if fd_thr > 0:
                        surf_func[hemi] = surf_func[hemi][idx_tp,:]
                    
                    # Regress out confounds (write a function for this)
                    b = np.linalg.lstsq(Z, surf_func[hemi], rcond=None)[0]
                    surf_func[hemi] = surf_func[hemi] - np.matmul(Z, b)
        
                    # For each cortical parcel, extract the average timecourse
                    U = np.unique(annot[hemi])
                    surf_parc[hemi] = np.zeros((surf_func[hemi].shape[0], U.shape[0]))
                    for parc in U:
                        surf_parc[hemi][:,parc] = np.mean(surf_func[hemi][:,annot[hemi] == parc], axis=1)
                     
                # Load the subcortical segmentation (FS "aseg" files) in MNI space
                [aseg, aseg_siz, aseg_aff]   = loadNifti(os.path.join(fmriprepdir, subj, ses, 'func', '{}{}_task-rest{}_space-T1w_desc-aseg_dseg.nii.gz'.format(subj, ses if not ses else '_{}'.format(ses), run if not run else '_{}'.format(run))))
                
                # Load their functional data in MNI space (volume)
                [vol_func, vol_siz, vol_aff] = loadNifti(os.path.join(fmriprepdir, subj, ses, 'func', '{}{}_task-rest{}_space-T1w_desc-preproc_bold.nii.gz'.format(subj, ses if not ses else '_{}'.format(ses), run if not run else '_{}'.format(run)))) 
                
                if remove_vols > 0:
                    vol_func = vol_func[remove_vols:,:]

                if fd_thr > 0:
                    vol_func = vol_func[idx_tp,:]

                # Regress out confounds
                b = np.linalg.lstsq(Z, vol_func, rcond=None)[0]
                vol_func = vol_func - np.matmul(Z, b)
                   
                # For each subcortical parcel, extract the average timecourse
                vol_parc = {}
                for hemi in ['L', 'R']:
                    vol_parc[hemi] = np.zeros((vol_func.shape[0], len(aseg_list[hemi])))
                    for pidx, parc in enumerate(aseg_list[hemi]):
                        vol_parc[hemi][:,pidx] = np.mean(vol_func[:,np.squeeze(aseg == parc)], axis=1)
                
                # Compute the "netmat" between the regions (cortical and subcortical),
                # and unwrap it
                all_parc = np.concatenate((surf_parc['L'][:,1:], surf_parc['R'][:,1:], vol_parc['L'], vol_parc['R']), axis=1)
                rmat = np.corrcoef(all_parc, rowvar=False);
                if partial == True:
                    rinv = np.linalg.pinv(rmat)
                    diag = np.diagonal(rinv)[:,None]
                    rmat = np.multiply(-rinv, np.power(np.multiply(diag, diag.T), -.5))
                uw = unWrap(rmat)
                uw = np.arctanh(uw)
                
                # Save the unwrapped as a text file (row only)
                if not os.path.exists(os.path.join(netmatsdir, subj, ses, 'func')):
                    os.makedirs(os.path.join(netmatsdir, subj, ses, 'func'))
                np.savetxt(os.path.join(netmatsdir, subj, ses, 'func', '{}{}_task-rest{}_netmat-{}_atlas-Schaefer2018-200P+17N_space-T1w.csv'.format(subj, ses if not ses else '_{}'.format(ses), run if not run else '_{}'.format(run), 'partial' if partial else 'full')), uw, delimiter=",")
                np.savetxt(os.path.join(netmatsdir, subj, ses, 'func', '{}{}_task-rest{}_netmat-{}_atlas-Schaefer2018-200P+17N_space-T1w.txt'.format(subj, ses if not ses else '_{}'.format(ses), run if not run else '_{}'.format(run), 'partial' if partial else 'full')), np.array([Z.shape[0]]).astype('int'), delimiter=",")
                imageio.imwrite(os.path.join(netmatsdir, subj, ses, 'func', '{}{}_task-rest{}_netmat-{}_atlas-Schaefer2018-200P+17N_space-T1w.png'.format(subj, ses if not ses else '_{}'.format(ses), run if not run else '_{}'.format(run), 'partial' if partial else 'full')), (255*(rmat+1)/2).astype('uint8'))
                plt.imshow(rmat)
                plt.show()
