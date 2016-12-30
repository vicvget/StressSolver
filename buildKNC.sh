rm -r CMakeFiles
rm -r CMakeCache.txt
cmake -D USE_KNC=true ""
make -j 8
