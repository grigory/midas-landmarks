#!/bin/bash
ulimit -s 512000
./mapgen map.co map.gr 1000000 > map.out
