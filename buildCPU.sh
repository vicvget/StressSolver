rm -r CMakeFiles
rm -r CMakeCache.txt
cmake -D USE_AVX=true ""
make -j 8

