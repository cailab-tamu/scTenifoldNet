
%% add everything to the path
path0 = pwd;
addpath(genpath(path0));

%% random seed generation
rand('twister', sum(100 * clock));

%% build all the mex files every time (dumb)
cd 'Alignment/dtw';
mex dtwFord.cpp;
mex dtwBack.cpp;
cd(path0);

cd 'support/';
mex rowBd.cpp;
mex -largeArrayDims knnsearch_same.cpp
cd(path0);

clear p paths path0 prefix;
