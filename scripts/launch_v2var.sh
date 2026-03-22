#!/bin/bash
NVRTC_PATH=$(python3 -c "import nvidia.cuda_nvrtc, os; print(os.path.dirname(nvidia.cuda_nvrtc.__file__)+'/lib')")
export LD_LIBRARY_PATH="${NVRTC_PATH}:${LD_LIBRARY_PATH}"
export CC=/usr/bin/gcc
nohup python3 /workspace/scripts/v2_variance_test.py > /workspace/v2var.log 2>&1 &
echo "launched pid=$!"
