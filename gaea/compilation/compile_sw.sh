#!/bin/bash
# Script to compile the non hydrostatic solver

if [ "$#" -ne 2 ]; then
  echo "Usage: ./compile.sh COMP clean"
  echo "choices for COMP :  debug, repro, prod"
  echo "choices for clean:  noclean, clean, cleanall"
  exit 1
fi

COMP="$1"
clean="$2"

model=sw
cd ../../SHiELD_build/Build
./COMPILE solo $model $COMP intel 64bit $clean

# Copy the log file to test directory
cp -r build_solo_$model.$COMP.64bit.intel.out ../../gaea/compilation/build_solo_$model.$COMP.64bit.intel.out
