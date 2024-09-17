#!/bin/bash
source /cms/cmsset_default.sh
echo $VO_CMS_SW_DIR
export SCRAM_ARCH=el9_amd64_gcc12
WORK=$workpath; shift
SRC=$workpath; shift
cd $SRC; 
echo "Getting ENV from $SRC"
eval $(scramv1 runtime -sh);
echo "Running in $WORK"
cd $WORK;
ulimit -c 0
echo "Will execute $command"
exec $command
