%clear all;
chdir('C:\OSU3\Repo\bundle-adjustment');
eigen_inc = ['-I' 'C:\OSU3\Repo\bundle-adjustment\Eigen'];

mex(eigen_inc, 'ba_algo.cpp', 'optim.cpp', 'matlab.cpp', 'core.cpp');
mex(eigen_inc, 'backproject.cpp', 'optim.cpp', 'matlab.cpp', 'core.cpp');
%mex(eigen_inc, 'multiray_triangulate.cpp', 'optim.cpp', 'matlab.cpp', 'core.cpp');


