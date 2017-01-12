export SINK_LD_LIBRARY_PATH=$MIC_LIBRARY_PATH
micnativeloadex _Debug/bin/test -a "64 1 0" -e "KMP_AFFINITY=verbose OMP_NUM_THREADS=10"
