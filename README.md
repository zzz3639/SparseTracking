# SparseTracking

This package is one of the supportive materials of paper ¡±SparseTracking: a super-high density particle tracking algorithm for single molecule localization microscopy¡±. It is designed for fluorescent emitter(molecule) tracking at super-high density. Other functions such as simulation of microscopic image series and localization of molecule motion barriers from trajectories are also included.

We write the core ¡±sparse linking EM¡± (SLEM) algorithm using C language, and we write other functions using matlab. The C codes need to be compiled before the use of this package, any compiler support C99 standard(gcc for linux, intel c/c++ compiler or mingw for windows) should be OK to do this compile.

Our package doesn¡¯t require extra software or package to execute except for ¡±SparseTracking Visualize¡±. ¡±SparseTracking Visualize¡± is written in the aim of visualizing trajectories, as well as reading/writing ¡±tif¡± files and movie files.

Find "manual.pdf" under folder "doc". 

We will upload more demos and we will continue to supplement "manual.pdf" in recent days.
