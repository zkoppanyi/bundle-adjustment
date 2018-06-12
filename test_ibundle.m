clear variables; clc; close all;

KNOWN = 1;
UNKNOWN = 2;

%% Settings

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

scale = 4.87e-6;

% Camera with distortion
% cam1 = [0.05 -0.1*1e-4 1e-6 3.5*1e-4 1*1e-6 0 -2.5*1e-5 -2.5*1e-5];
% cam2 = [0.02  0.1*1e-4 1e-6 3.5*1e-4 1.3*1e-6 0 -2.5*1e-5 -2.5*1e-5];
% cams = [2 2 cam1 KNOWN; 
%         3 2 cam2 KNOWN];

% Camera without distortion
% cam1 = [0.05 -0.1*1e-4 1e-6];
% cam2 = [0.02  0.1*1e-4 1e-6];
cam1 = [0.05  0 0];
cam2 = [0.05  0 0];
cams = [2 1 cam1 KNOWN; 
        3 1 cam2 KNOWN];

% generate tie points as grid
objx = -25 : 5 : 25;
objy = -25 : 5 : 25;
[objx, objy] = meshgrid(objx, objy);
objx = objx(:);
objy = objy(:);
n_tie_pts = length(objx);

% generate control points as grid
% obj_ct_x = -25 : 15 : 25;
% obj_ct_y = -25 : 15 : 25;
% [obj_ct_x, obj_ct_y] = meshgrid(obj_ct_x, obj_ct_y);
% obj_ct_x = obj_ct_x(:);
% obj_ct_y = obj_ct_y(:);
% objx = [objx; obj_ct_x];
% objy = [objy; obj_ct_y];
% n_ct_pts = length(obj_ct_x);

n_test_pt = length(objx);
obj_pts = [(1:n_test_pt)', objx, objy, (rand(n_test_pt, 1)-0.5)*10, ones(n_test_pt, 1)*KNOWN];

imgs = [1 0  0  56     -3/180*pi  1/180*pi 85/180*pi 2 KNOWN;
        2 50  50  56   -3/180*pi -5/180*pi 95/180*pi 2 KNOWN;
%         3 0  60 70   -3/180*pi 5/180*pi 75/180*pi 2 KNOWN;
%         4 65 55 50   30/180*pi 5/180*pi 55/180*pi 2 KNOWN;
%         5 50 20 87  -10/180*pi 5/180*pi 55/180*pi 3 KNOWN;
%         6 20 50 56   -3/180*pi 5/180*pi 55/180*pi 3 KNOWN;
%         7 50 50 70   -3/180*pi 5/180*pi 75/180*pi 3 KNOWN;
%         8 55 55 20   30/180*pi 5/180*pi 55/180*pi 3 KNOWN
];

    
%% Initialization with backprojection
disp('Starting backproject...');
tic

%imgs = [1 23 15 100  -3/180*pi 5/180*pi 95/180*pi 1 2];
img_pts = [1 0 0 1 1 1];
img_pts = backproject(img_pts, imgs, obj_pts, cams);
cams(:,2) = 1;

%plot_problem(1, img_pts, imgs, obj_pts, cams, 'r');

toc
disp('End backproject.')

%% Incremental bundle adjustment
img_pts1 = img_pts(img_pts(:,4) == 1, :);
img_pts2 = img_pts(img_pts(:,4) == 2, :);
img_pts1(:,2:3) = img_pts1(:,2:3) + [600 400];
img_pts2(:,2:3) = img_pts2(:,2:3) + [600 400];
F = estimateFundamentalMatrix(img_pts1(:,2:3), img_pts2(:,2:3), 'Method', 'Norm8Point');

K1 = cameraIntrinsics(cam1(1),cam1(2:3) + [600 400], [600 400]);
K2 = cameraIntrinsics(cam2(1),cam2(2:3) + [600 400], [600 400]);
[rot, t] = relativeCameraPose(F,K1,K2, img_pts1(:,2:3), img_pts2(:,2:3));
eul_ang = rotm2eul(rot);
eul_ang/pi*180
t*sqrt(50^2+50^2)

%stereoParams = stereoParameters(K1,K2,rot,t)

%img_pts(:,2:3) = img_pts(:,2:3) + normrnd(0, 0.0001, size(img_pts, 1), 2);

obj_pts0 = obj_pts;
tie_pts_idx  = 1 : n_tie_pts;
%obj_pts0(tie_pts_idx, 2:5) = repmat([0 0 0 UNKNOWN], length(tie_pts_idx), 1);
obj_pts0(tie_pts_idx, end) = UNKNOWN;
%obj_pts0(tie_pts_idx, 2:4) = obj_pts(tie_pts_idx, 2:4)+3;
control_pts_idx = find(obj_pts0(:, end) == KNOWN);

fprintf("Number of tie points: %i\n", length(tie_pts_idx));
fprintf("Number of control points: %i\n", length(control_pts_idx));


cam10 = [0.05 0 0];
cam20 = [0.02 0 0];
cams0 = [2 2 cam10 UNKNOWN; 
         3 2 cam20 UNKNOWN];
   
imgs0 = imgs;
%imgs0(:, 2:4) = imgs0(:, 2:4) + rand(size(imgs0,1), 3)*5;
%imgs0(:, 5:7) = imgs0(:, 5:7) + rand(size(imgs0,1), 3)*5/180*pi;
imgs0(1, 2:4) = 0;
imgs0(2, 5:7) = 0;
imgs0(1, 2:4) = [0 0 56];
imgs0(2, 2:7) = [t(1:2)*sqrt(50^2+50^2) 56 eul_ang(3) eul_ang(2) eul_ang(1)];
imgs0(:, end) = UNKNOWN;

disp('Starting bundle...');
tic

% 1. First solve without distortions
cams0(:,2) = 1;
sol = ba_algo(img_pts, imgs0, obj_pts0, cams0);
cams0(:, 3:5) = sol.cams(:, 3:5);
%sol.cams(:,2) = 1;

%% Visualization
plot_problem(1, sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, 'g');


