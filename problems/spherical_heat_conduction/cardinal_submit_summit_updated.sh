#/bin/bash

# =============================================================================
# Cardinal-specific
# =============================================================================

: ${PROJ_ID:=""}
: ${NEKRS_HOME:="$(realpath ../..)"}
: ${OCCA_CACHE_DIR:="$PWD/.cache/occa"}
: ${CARDINAL_BIN:=../../cardinal-opt}
: ${CARDINAL_APP:=nek}
: ${CARDINAL_I:=nek.i}

export OPENMC_CROSS_SECTIONS=/gpfs/alpine/nfi114/proj-shared/rahaman/cross_sections/nndc_hdf5/cross_sections.xml
export OPENMC_ENDF_DATA=/gpfs/alpine/nfi114/proj-shared/rahaman/cross_sections/endf-b-vii.1

# =============================================================================
# Adapted from nrsqsub_summit
# =============================================================================

NVME_HOME="/mnt/bb/$USER/"
XL_HOME="/sw/summit/xl/16.1.1-3/xlC/16.1.1"

: ${CPUONLY:=0}
export NEKRS_HOME
export OCCA_CACHE_DIR
export NEKRS_HYPRE_NUM_THREADS=1
export OGS_MPI_SUPPORT=1
export OCCA_CXX="$XL_HOME/bin/xlc" 
export OCCA_CXXFLAGS="-O3 -qarch=pwr9 -qhot -DUSE_OCCA_MEM_BYTE_ALIGN=64" 
export OCCA_LDFLAGS="$XL_HOME/lib/libibmc++.a"

#export OCCA_VERBOSE=1
#export OMPI_LD_PRELOAD_POSTPEND=$OLCF_SPECTRUM_MPI_ROOT/lib/libmpitrace.so

#export PAMI_ENABLE_STRIPING=1
#export PAMI_IBV_ADAPTER_AFFINITY=1
#export PAMI_IBV_DEVICE_NAME="mlx5_0:1,mlx5_3:1"
#export PAMI_IBV_DEVICE_NAME_1="mlx5_3:1,mlx5_0:1"

if [ -z "$PROJ_ID" ]; then
  echo "ERROR: PROJ_ID is empty"
  exit 1
fi

if [ $# -ne 3 ]; then
  echo "usage: [PROJ_ID] [CPUONLY=1] $0 <casename> <number of compute nodes> <hh:mm>"
  exit 0
fi

module unload darshan-runtime
module load gcc

bin=$CARDINAL_BIN
case=$1
nodes=$2
gpu_per_node=6
cores_per_socket=21
let nn=$nodes*$gpu_per_node
let ntasks=nn
time=$3
backend=CUDA

if [ $CPUONLY -eq 1 ]; then
  backend=CPU
  let nn=2*$nodes
  let ntasks=$nn*$cores_per_socket
fi 

if [ ! -f $bin ]; then
  echo "Cannot find" $bin
  exit 1
fi

if [ ! -f $case.par ]; then
  echo "Cannot find" $case.par
  exit 1
fi

if [ ! -f $case.co2 ]; then
  echo "Cannot find" $case.co2
  exit 1
fi

if [ ! -f $case.udf ]; then
  echo "Cannot find" $case.udf
  exit 1
fi

if [ ! -f $case.oudf ]; then
  echo "Cannot find" $case.oudf
  exit 1
fi

if [ ! -f $case.re2 ]; then
  echo "Cannot find" $case.re2
  exit 1
fi

mkdir -p $OCCA_CACHE_DIR 2>/dev/null

while true; do
  read -p "Do you want precompile (recommended)? [N]" yn
  case $yn in
    [Yy]* )
      echo $NEKRS_HOME
      mpirun -pami_noib -np 1 $NEKRS_HOME/bin/nekrs --setup $case --build-only $ntasks --backend $backend;
      if [ $? -ne 0 ]; then
        exit 1
      fi
      break ;;
    * )
      break ;;
  esac
done

if [ $CPUONLY -eq 1 ]; then
  jsrun="jsrun -X 1 -n$nodes -r1 -a1 -c1 -g0 -b packed:1 -d packed cp -a $OCCA_CACHE_DIR/* $NVME_HOME; export OCCA_CACHE_DIR=$NVME_HOME; jsrun -X 1 -n$nn -a$cores_per_socket -c$cores_per_socket -g0 -b packed:1 -d packed $bin --app $CARDINAL_APP -i $CARDINAL_I --nekrs-setup $case"
else
  jsrun="jsrun -X 1 -n$nodes -r1 -a1 -c1 -g0 -b packed:1 -d packed cp -a $OCCA_CACHE_DIR/* $NVME_HOME; export OCCA_CACHE_DIR=$NVME_HOME; jsrun --smpiargs='-gpu' -X 1 -n$nn -r$gpu_per_node -a1 -c2 -g1 -b rs -d packed $bin --app $CARDINAL_APP -i $CARDINAL_I --nekrs-setup $case" 
fi

cmd="bsub -nnodes $nodes -alloc_flags NVME -W $time -P $PROJ_ID -J cardinal_$case \"${jsrun}\""
echo $cmd
$cmd
