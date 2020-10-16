# Linke_jointCCAICA
Code used in Linke et al., Biol Psychiatry

These series of scripts assumes that the data is named and organized according to the Brain Imaging Data Structure (BIDS). For details see: https://bids-specification.readthedocs.io/en/stable/ 

Before we start, we create a virtual environment by executing the a01_create_virtualenv.sh script. This needs to be done only once. 
./a01_create_virtualenv.sh

Thereafter before running any of the scripts the environment should be activated:
source source.this

Next, we run MRIQC by exectuing a02_mriqc. This script assumes that in the IDs of the participants that should be processed live in $ROOTDIR/listst/MYLIST_NIMH.txt and MYLIST_CMI.txt

Based on the MRIQC output images that contain nonsteady state are removed from the timeseries in the BIDS directories. The respective numbers need to be entered into a03_remove_nonsteadystate.sh (line 3). DO NOT RUN TWICE!!!

Now we run a04_fmriprep.sh, which will preprocess the functional and structural data as described in Linke et al., 2020.

We are ready to create the connectivity matrices using b01_netmats.py (call spyder and execute). As the some of the CMI data had be acquired in two runs, we need to merge them by applying b02_merge_runs (call matlab and execute there).

The matrices and the clinical data can now be analyzed with c01_ccaica.m.
