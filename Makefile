all: conversion.pdf conversion-full.pdf

conversion.pdf: conversion.tex
	pdflatex conversion.tex

conversion-full.pdf: conversion-full.tex
	pdflatex conversion-full.tex
	bibtex conversion-full
	pdflatex conversion-full.tex
	pdflatex conversion-full.tex
	cp conversion-full.pdf "/home/salvi/pCloudDrive/Public Folder/"
