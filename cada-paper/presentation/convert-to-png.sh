#!/bin/bash

convert -density 600 slides.pdf slides-%02d.png
cp slides-00.png slides-22.png # dummy slide at the end
