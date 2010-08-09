#!/bin/bash
make
./midas.exe -generate $1 -pool 20000 -randc 10 -upper 10 -points 2048 -t1 800 -t2 990 -radius 1000 -seed 123 > landmark.lmk
./midas.exe -evaluate $1 -seed 566 < landmark.lmk
#./midas.exe -evaluate bay.gr -seed 346 < landmark.lmk
#./midas.exe -evaluate bay.gr -seed 610 < landmark.lmk
#./midas.exe -evaluate bay.gr -seed 22 < landmark.lmk

