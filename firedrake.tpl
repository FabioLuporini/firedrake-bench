#!/bin/bash --login

#PBS -N %(jobname)s
#PBS -l walltime=%(walltime)s:0
#PBS -l select=%(nodes)d
#PBS -A n02-NEK006789
#PBS -q standard
#PBS -m eba
#PBS -M david.ham@imperial.ac.uk

#set -v
set -x

# Make sure any symbolic links are resolved to absolute path 
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

echo Running in $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

LOGFILE=%(jobname)s.${PBS_JOBID}.log

module swap PrgEnv-cray PrgEnv-gnu
module add anaconda
module use /work/n02/n02/dham/fdrake/modules
module load firedrake

export FIREDRAKE_FFC_KERNEL_CACHE_DIR=$WORK/firedrake-cache
#export FIREDRAKE_FFC_KERNEL_CACHE_DIR=$WORK/firedrake-cache/${PBS_JOBID}
#export FIREDRAKE_FFC_KERNEL_CACHE_DIR=$TMPDIR
export PYOP2_LAZY=0
export PYOP2_PROFILING=1
export PYOP2_BACKEND_COMPILER=gnu
export PYOP2_SIMD_ISA=avx
export PYOP2_CACHE_DIR=$WORK/pyop2-cache
#export PYOP2_CACHE_DIR=$WORK/pyop2-cache/${PBS_JOBID}
#export PYOP2_CACHE_DIR=$TMPDIR
#export LD_LIBRARY_PATH=$ANACONDA_LIB:$LD_LIBRARY_PATH
#export PYTHONPATH=$WORK/firedrake-bench:$WORK/pybench:$WORK/firedrake:$WORK/PyOP2:$PYTHONPATH
export PYTHONPATH=$WORK/firedrake-bench:$WORK/pybench:$PYTHONPATH
export PETSC_OPTIONS=-log_summary
# Prevent matplotlib from accessing /home
export HOME=$WORK
export XDG_CONFIG_HOME=''

# MPI (man intro_mpi)
export MPICH_NEMESIS_ASYNC_PROGRESS=MC
export MPICH_MAX_THREAD_SAFETY=multiple
export MPICH_CPUMASK_DISPLAY=1

module load firedrake

echo -n Started at | tee -a $LOGFILE
date | tee -a $LOGFILE

echo | tee -a $LOGFILE
echo Running %(script)s %(args)s 2>&1  | tee -a $LOGFILE
echo | tee -a $LOGFILE
echo aprun -b -n %(ptotal)d -N %(ppn)d -S %(pnuma)d python %(script)s.py 2>&1  | tee -a $LOGFILE
#echo aprun -b -n %(ptotal)d -N %(ppn)d -S %(pnuma)d python %(script)s.py %(args)s 2>&1  | tee -a $LOGFILE


python %(script)s.py %(args)s 2>&1  | tee -a $LOGFILE
aprun -b -n %(ptotal)d -N %(ppn)d -S %(pnuma)d python %(script)s.py %(args)s 2>&1  | tee -a $LOGFILE

echo -n Finished at | tee -a $LOGFILE
date | tee -a $LOGFILE
