#!/bin/bash

SWARMFILE=${ROOTDIR}/slurm/swarm_fmriprep.txt
echo -n "" > ${SWARMFILE}
for DATASET in _CMI _NIMH; do  
for s in $(cat ${ROOTDIR}/lists/list${DATASET}.csv) ; do
   if [[ ! -e ${ROOTDIR}/derivatives${DATASET}/fmriprep/sub-${s}.html ]] || \
      [[ $(grep "No errors to report" ${ROOTDIR}/derivatives${DATASET}/fmriprep/sub-${s}.html |wc -l ) -ne 1 ]] ; then 
      echo "let \"rnd = ${RANDOM} % 300\" ; sleep \${rnd} ; \
            export TMPDIR=/lscratch/\${SLURM_JOBID} ; \
            mkdir -p \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk ; \
            fmriprep ${ROOTDIR}/BIDS${DATASET} \${TMPDIR}/${s}.out participant --participant_label ${s} -w \${TMPDIR}/${s}.wrk --use-aroma --notrack --output-space T1w template fsnative fsaverage --nthreads 1 --omp-nthreads 1 --skip_bids_validation --notrack ; \
            mkdir -p ${ROOTDIR}/derivatives${DATASET}+/ ; \
            rsync -a \${TMPDIR}/${s}.out/ ${ROOTDIR}/derivatives${DATASET}+/ ; \
            rm -rf \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk " >> ${SWARMFILE}
   fi
done
done
swarm -f ${SWARMFILE} -g 10 -t 1 -p 1 --gres=lscratch:10 --logdir ${ROOTDIR}/slurm --time=120:00:00 --job-name merfp

