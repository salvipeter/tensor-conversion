#!/bin/bash

ffmpeg -i slides.mp4 -i slides.mp3 -codec copy -shortest presentation.mp4
