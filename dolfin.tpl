#!/bin/bash --login

#PBS -N %(jobname)s
#PBS -l walltime=%(walltime)s:0
#PBS -l select=%(nodes)d
#PBS -A n02-NEK006789
#PBS -m eba
#PBS -M david.ham@imperial.ac.uk

#set -v
set -x

# Make sure any symbolic links are resolved to absolute path 
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

echo Running in $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

LOGFILE=%(jobname)s.${PBS_JOBID}.log

module load dot
module use /work/y07/y07/fdrake/modules
#module load fenics/master

export INSTANT_CACHE_DIR=$WORK/.instant/cache
export INSTANT_ERROR_DIR=$WORK/.instant/error
export TMPDIR=$WORK/.instant

export PYTHONPATH=$WORK/pybench:$PYTHONPATH
export PETSC_OPTIONS=-log_summary
# Prevent matplotlib from accessing /home
export HOME=$WORK
# MPI (man intro_mpi)
export MPICH_NEMESIS_ASYNC_PROGRESS=MC
export MPICH_MAX_THREAD_SAFETY=multiple
export MPICH_CPUMASK_DISPLAY=1

module load fenics
module list
env
echo -n Started at | tee -a $LOGFILE
date | tee -a $LOGFILE

echo | tee -a $LOGFILE
echo Running %(script)s %(args)s 2>&1  | tee -a $LOGFILE
echo | tee -a $LOGFILE

aprun -n %(ptotal)d -N %(ppn)d -S %(pnuma)d python %(script)s.py %(args)s 2>&1  | tee -a $LOGFILE

echo -n Finished at | tee -a $LOGFILE
date | tee -a $LOGFILE
