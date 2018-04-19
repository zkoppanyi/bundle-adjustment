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
%cam1 = [0.05 -0.1*1e-4 1e-6 3.5*1e-8 1*1e-12 0 -2.5*1e-8 -2.5*1e-8];
%cam2 = [0.02  0.1*1e-4 1e-6 3.5*1e-8 1.3*1e-12 0 -2.5*1e-8 -2.5*1e-8];
cam1 = [0.05 -0.1*1e-4 1e-6 3.5*1e-4 1*1e-6 0 -2.5*1e-5 -2.5*1e-5];
cam2 = [0.02  0.1*1e-4 1e-6 3.5*1e-4 1.3*1e-6 0 -2.5*1e-5 -2.5*1e-5];
cams = [2 2 cam1 KNOWN; 
        3 2 cam2 KNOWN];

% Camera without distortion
% cam1 = [0.05 -0.1*1e-4 1e-6];
% cam2 = [0.02  0.1*1e-4 1e-6];
% cams = [2 1 cam1 KNOWN; 
%         3 1 cam2 KNOWN];

% generate tie points as grid
objx = -25 : 5 : 25;
objy = -25 : 5 : 25;
[objx, objy] = meshgrid(objx, objy);
objx = objx(:);
objy = objy(:);
n_tie_pts = length(objx);

% generate control points as grid
obj_ct_x = -25 : 15 : 25;
obj_ct_y = -25 : 15 : 25;
[obj_ct_x, obj_ct_y] = meshgrid(obj_ct_x, obj_ct_y);
obj_ct_x = obj_ct_x(:);
obj_ct_y = obj_ct_y(:);
objx = [objx; obj_ct_x];
objy = [objy; obj_ct_y];
n_ct_pts = length(obj_ct_x);

n_test_pt = length(objx);
obj_pts = [(1:n_test_pt)', objx, objy, (rand(n_test_pt, 1)-0.5)*10, ones(n_test_pt, 1)*KNOWN];

imgs = [1 0  0  87   -3/180*pi 5/180*pi 95/180*pi 2 KNOWN;
        2 0  0  56   -3/180*pi 5/180*pi 95/180*pi 2 KNOWN;
        3 0  60 70   -3/180*pi 5/180*pi 75/180*pi 2 KNOWN;
        4 65 55 50   30/180*pi 5/180*pi 55/180*pi 2 KNOWN;
        5 50 20 87  -10/180*pi 5/180*pi 55/180*pi 3 KNOWN;
        6 20 50 56   -3/180*pi 5/180*pi 55/180*pi 3 KNOWN;
        7 50 50 70   -3/180*pi 5/180*pi 75/180*pi 3 KNOWN;
        8 55 55 20   30/180*pi 5/180*pi 55/180*pi 3 KNOWN];

    
%% Initialization with backrpojection
disp('Starting backproject...');
tic

%imgs = [1 23 15 100  -3/180*pi 5/180*pi 95/180*pi 1 2];
img_pts = [1 0 0 1 1 1];
img_pts = backproject(img_pts, imgs, obj_pts, cams);
cams(:,2) = 1;

img_pts_ud = backproject(img_pts, imgs, obj_pts, cams);
cams(:,2) = 2;

toc
disp('End backproject.')

%% Bundle adjustment

%img_pts(:,2:3) = img_pts(:,2:3) + normrnd(0, 0.0001, size(img_pts, 1), 2);

obj_pts0 = obj_pts;
tie_pts_idx  = 1 : n_tie_pts;
%obj_pts0(tie_pts_idx, 2:5) = repmat([0 0 0 UNKNOWN], length(tie_pts_idx), 1);
obj_pts0(tie_pts_idx, end) = UNKNOWN;
obj_pts0(tie_pts_idx, 2:4) = obj_pts(tie_pts_idx, 2:4)+3;
control_pts_idx = find(obj_pts0(:, end) == KNOWN);

fprintf("Number of tie points: %i\n", length(tie_pts_idx));
fprintf("Number of control points: %i\n", length(control_pts_idx));


cam10 = [0.052 0 0 1e-4 1e-6 0 1e-5 1e-5];
cam20 = [0.012 0 0 1e-4 1e-6 0 1e-5 1e-5];
cams0 = [2 2 cam10 UNKNOWN; 
         3 2 cam20 UNKNOWN];
   
imgs0 = imgs;
imgs0(:, 2:4) = imgs0(:, 2:4) + rand(size(imgs0,1), 3)*5;
imgs0(:, 5:7) = imgs0(:, 5:7) + rand(size(imgs0,1), 3)*5/180*pi;
imgs0(:, end) = UNKNOWN;

disp('Starting bundle...');
tic

% 1. First solve without distortions
% cams0(:,2) = 1;
% sol = ba_algo(img_pts, imgs0, obj_pts0, cams0);
% cams0(:, 3:5) = sol.cams(:, 3:5);
% %sol.cams(:,2) = 1;


% 2. Solve with distortions
cams0(:,2) = 2;
%[sol, stoch] = ba_algo(img_pts, sol.imgs, obj_pts0, cams0);
%sol = ba_algo(img_pts, sol.imgs, obj_pts0, cams0);
sol = ba_algo(img_pts, imgs0, obj_pts0, cams0);

img_pts_ud2 = backproject(img_pts, sol.imgs, sol.obj_pts, sol.cams);
plot_problem(1, sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, 'g');

toc
disp('End bundle.')

%% Visualization

% Radial distrotion
idx = find(img_pts(:, 4) == 1);
dr = img_pts(idx,2:3) - img_pts_ud2(idx,2:3);
dr = sqrt(dr(:,1).^2 + dr(:,2).^2);
r = sqrt(img_pts_ud2(idx,2).^2 + img_pts_ud2(idx,3).^2);
figure(2); clf; hold on;
plot(r, dr, 'r*');
title('Radial distrotion');

show_radial_distortion_profile(sol, 0);

% Image plane
figure(1); clf; hold on;
idx = find(img_pts(:, 4) == 1);
plot(img_pts(idx,2), img_pts(idx,3), 'r.', 'MarkerSize', 10);
plot(img_pts_ud(idx,2), img_pts_ud(idx,3), 'b.', 'MarkerSize', 10);
plot(img_pts_ud2(idx,2), img_pts_ud2(idx,3), 'go');
legend('Distorted', 'Undistorted (ground truth)', 'Undistorted (solution)');
title('Image plane');
grid on;
axis equal;

% Stochastics results
% diagMxx = diag(stoch.Mxx);
% coorMxx = sqrt(diagMxx(1:3));
% rotMxx = sqrt(diagMxx(4:6));
% apostMll = sqrt(diag(stoch.Mll));
% fprintf('Coordinates apost. error: %.3f m\n', norm(coorMxx) );
% fprintf('  Rotations apost. error: %.3f deg\n', norm(rotMxx)/pi*180 );
% fprintf(' Img. coor. apost. error: %.3f pixel\n', norm(apostMll));

%
figure(3);
plot_problem(1, sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, 'g');

% Sparsity pattern
% figure(2);
% show_sparisty_pattern(stoch.J, 10);

%%
% J = stoch.J;
% H = stoch.J'*stoch.J;
% nH = size(H,1);
% n1 = nH - n_tie_pts*3+1;
% C = H(n1:nH, n1:nH);
% invC = inv(C);
% 
% nimgs = size(imgs,1)*3 + size(cams,1)*8*2;
% B = H(1:nimgs, 1:nimgs);
% 
% E1 = H(1:nimgs, n1:nH);
% E2 = H(n1:nH, 1:nimgs);
% 
% sum(sum(E1-E2')) == 0; %!!!
% E = E1;
% 
% S = B - E*invC*E';
% show_sparisty_pattern(S, 10);


