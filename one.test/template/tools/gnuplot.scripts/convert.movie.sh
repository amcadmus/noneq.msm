#!/bin/bash

convert fig-2d-tog.gif test%05d.png
ffmpeg -r 4 -i test%05d.png -y -an fig-2d-tog.mp4
rm -f test*.png
