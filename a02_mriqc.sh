#!/bin/bash

SWARMFILE=${ROOTDIR}/slurm/swarm_mriqc.txt
echo -n "" > ${SWARMFILE}
for DATASET in _NIMH _CMI ; do 
for s in $(cat ${ROOTDIR}/lists/MYLIST${DATASET}.txt) ; do
   if [[ ! -d ${ROOTDIR}/derivatives${DATASET}/mriqc/sub-${s} ]] ; then 
      echo "let \"rnd = ${RANDOM} % 300\" ; sleep \${rnd} ; \
            export TMPDIR=/lscratch/\${SLURM_JOBID} ; \
            mkdir -p \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk ; \
            mriqc ${ROOTDIR}/BIDS${DATASET} \${TMPDIR}/${s}.out participant --participant_label ${s} -w \${TMPDIR}/${s}.wrk --n_procs 1 --no-sub ; \
            mkdir -p ${ROOTDIR}/derivatives${DATASET}+/mriqc ; \
            rsync -a \${TMPDIR}/${s}.out ${ROOTDIR}/derivatives${DATASET}+/mriqc ; \
            rm -rf \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk " >> ${SWARMFILE}
   fi
done
done
swarm -f ${SWARMFILE} -g 40 -t auto --gres=lscratch:40 --logdir ${ROOTDIR}/slurm --time=120:00:00 --job-name mriqc

