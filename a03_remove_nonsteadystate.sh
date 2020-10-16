#!/bin/bash

BIDS_LIST=(BIDS_CMI BIDS_NIMH)
DUMMY_LIST=(5 4)

for ((i=0; i<=${#BIDS_LIST[@]}; i++)) ; do
bids=${BIDS_LIST[i]}
dummy_trs=${DUMMY_LIST[i]}
cd ${ROOTDIR}/${bids}/
for sub in sub-* ; do
   echo "=== [ ${sub} ] ==="
   cd ${ROOTDIR}/${bids}/${sub}
#   for ses in ses-* ; do
ses=""
     cd ${ROOTDIR}/${bids}/${sub}/${ses}/func
     for f in *_bold.nii.gz ; do
        echo ${f}
        if [[ -e ${f} ]] ; then
           dim4=$(${FSLDIR}/bin/fslinfo ${f} | grep ^dim4 | awk '{print $2}')
           newdim4=$(expr ${dim4} - ${dummy_trs})
           ${FSLDIR}/bin/fslroi ${f} ${f} ${dummy_trs} ${newdim4}
        fi
     done
#   done
done
done
