
HOSTNAME=$(uname -n | awk -F. '{print $1}')

#export PS1="\[\033[1;34m\]\u@\h:\w $\[\033[0m\] "
PROJECT_NAME=MErest
export PS1="[${PROJECT_NAME}] ${PS1}"
#if   [[ "${HOSTNAME}" == "felix" ]] || \
#     [[ "${HOSTNAME}" == "biowulf" ]] ; then
  ROOTDIR=/data/EDB/MErest
  module load afni
  module load dcm2niix
  module load fsl/5.0.11
  module load freesurfer/6.0.0
  module load R/3.5.2
  ANTSDIR=/data/EDB/opt/ANTs/2.3.1/bin
  ANTSPATH=/data/EDB/opt/ANTs/2.3.1/bin
  C3DDIR=/data/EDB/opt/c3d/1.0.0/bin
  AROMADIR=${ROOTDIR}/code/ICA-AROMA-0.4.4-beta
  CRN_SHARED_DATA=${ROOTDIR}/templates
  export PATH=${PATH}:${ANTSDIR}:${ANTSDIR}/../Scripts:${C3DDIR}:${AROMADIR}
  export export TEMPLATEFLOW_HOME=${ROOTDIR}/templateflow
  export ANTSDIR
  export ANTSPATH
  export CRN_SHARED_DATA
  source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
  source ${ROOTDIR}/code/env-hpc/bin/activate
  module load matlab
#elif [[ "${HOSTNAME}" == "springsteen" ]] || \
#     [[ "${HOSTNAME}" == "bobdylan" ]] || \
#     [[ "${HOSTNAME}" == "jsbach" ]] ; then
#  echo "Error: Not configured to run on ${HOSTNAME}."
#fi

export ROOTDIR
export PATH=${PATH}:${ROOTDIR}/code:${AROMADIR}
alias ls="ls --color=auto"
umask 0002

