%%
% Goal: Test backprojection code
%%

clear all; clc;

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

% Cameras
% ID, cam_type, f, cx, cy, type
cams = [1 0 1.0 0.0 0.0 1];

% Object points
% ID, X, Y, Z, type
n_test_pt = 10;
rand_pts = [(rand(n_test_pt, 1)-0.5)*150, (rand(n_test_pt, 1)-0.5)*150, (rand(n_test_pt, 1)-0.5)*30];
obj_pts = [(1:n_test_pt)' rand_pts ones(n_test_pt, 1)];

%Images 
% ID, X, Y, Z, omega, phi kappa, cam_id, type
imgs = [1 0 0 100 3/180*pi 5/180*pi 25/180*pi 1 1];

% Image points
% ID, x, y, img_id, obj_id, type
img_pts =[1 0 0 1 1 1];

img_pts = backproject(img_pts, imgs, obj_pts, cams);

opts = plotset;
opts.is_show_rays = 1;
plot_problem(1, img_pts, imgs, obj_pts, cams, opts);

