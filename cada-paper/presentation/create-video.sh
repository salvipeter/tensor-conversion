#!/bin/bash

ffmpeg -ts_from_file 2 -i 'slides-%02d.png' slides.mp4
