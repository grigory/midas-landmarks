#!/bin/bash
g++ mapgen.cpp -o mapgen
ulimit -s 512000
./mapgen map.co map.gr 1000000 > map.out
