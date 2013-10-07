#!/bin/bash

make -C tools/angles/ clean
make -j -C tools/angles/

./tools/angles/average.traj.metastable
