# This needs to be done just once
module load python/3.6

ROOTDIR=/MYPATH

virtualenv env-hpc
source ${ROOTDIR}/code/env-hpc/bin/activate
pip install sentry-sdk
pip install requests
pip install duecredit
pip install mriqc
pip install triangle
pip install fmriprep

