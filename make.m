%clear all;

% Compiler: Visual C++ 2015
mex COMPFLAGS="$COMPFLAGS /O2 /arch:AVX /openmp" -IC:\OSU3\Repo\bundle-adjustment\Eigen 'ba_problem.cpp' 'core.cpp' 'optim.cpp' 'matlab.cpp' 
mex COMPFLAGS="$COMPFLAGS /O2 /arch:AVX /openmp" -IC:\OSU3\Repo\bundle-adjustment\Eigen 'ba_algo.cpp' 'core.cpp' 'optim.cpp' 'matlab.cpp' 
%mex COMPFLAGS="$COMPFLAGS /O2 /arch:AVX /openmp" -IC:\OSU3\Repo\bundle-adjustment\Eigen 'backproject.cpp' 'core.cpp' 'optim.cpp' 'matlab.cpp' 
%mex COMPFLAGS="$COMPFLAGS" -IC:\OSU3\Repo\bundle-adjustment\Eigen 'core.cpp' 'optim.cpp' 'matlab.cpp' 'multiray_triangulate.cpp'


