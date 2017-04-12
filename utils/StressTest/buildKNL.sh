rm -r CMakeFiles
rm -r CMakeCache.txt
cmake -D USE_KNL=true ""
make -j 8
