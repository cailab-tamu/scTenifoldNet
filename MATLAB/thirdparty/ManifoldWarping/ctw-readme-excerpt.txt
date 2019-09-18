(Much of the support code for Manifold Warping is based on Feng Zhou's CTW implementation.
The original CTW readme file is excerpted below.)

---

This page contains software and instructions for Canonical Time Warping (CTW) [1]. All the functions have been written and documented in Matlab format. We additionally provide C++ implementations of some dynamic programming routines which involve many loops and are notoriously slow in Matlab.
    [1] F. Zhou and F. de la Torre, "Canonical time warping for alignment of human behavior", Neural Information Processing Systems (NIPS), 2009.

For each C++ file, we provide its corresponding Matlab version. For instance, you can use "dtwFordSlow.m" instead of "dtwFord.cpp". They have the same interface in both input and output.


This software is free for use in research projects. If you publish results obtained using this software, please use this citation:

@inproceedings{Zhou_2009_6478,
   author = {Feng Zhou and Fernando de la Torre},
   title = {Canonical Time Warping for Alignment of Human Behavior},
   booktitle = {Neural Information Processing Systems Conference (NIPS)},
   month = {December},
   year = {2009},
}

If you have any questions, please feel free to contact Feng Zhou (zhfe99@gmail.com).
