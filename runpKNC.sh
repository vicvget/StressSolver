export SINK_LD_LIBRARY_PATH=$MIC_LIBRARY_PATH
micnativeloadex _Debug/bin/test -a "64 10 0" -e "KMP_AFFINITY=verbose OMP_NUM_THREADS=118"
