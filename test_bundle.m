clear variables; clc; close all;

%% Settings

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

scale = 4.87e-6;
cam1 = [0.05 -0.1*1e-4 1e-6 3.5*1e-8 1*1e-12 0 -2.5*1e-8 -2.5*1e-8];
cam2 = [0.02  0.1*1e-4 1e-6 3.5*1e-8 1.3*1e-12 0 -2.5*1e-8 -2.5*1e-8];
cams = [2 2 cam1 1; 3 2 cam2 1];

% n_test_pt = 120;
% rand_pts = [(rand(n_test_pt, 1))*25, (rand(n_test_pt, 1))*25, (rand(n_test_pt, 1)-0.5)*30];
% obj_pts = [(1:n_test_pt)' rand_pts ones(n_test_pt, 1)];

objx = -25 : 5 : 25;
objy = -25 : 5 : 25;
[objx, objy] = meshgrid(objx, objy);
objx = objx(:);
objy = objy(:);
n_test_pt = length(objx);
obj_pts = [(1:n_test_pt)', objx, objy, (rand(n_test_pt, 1)-0.5)*10, ones(n_test_pt, 1)];

imgs = [1 0  0  87   -3/180*pi 5/180*pi 95/180*pi 2 1;
        2 0  0  56   -3/180*pi 5/180*pi 95/180*pi 2 1;
        3 0  60 70   -3/180*pi 5/180*pi 75/180*pi 2 1;
        4 65 55 50   30/180*pi 5/180*pi 55/180*pi 2 1;
        5 50 20 87  -10/180*pi 5/180*pi 55/180*pi 3 1;
        6 20 50 56   -3/180*pi 5/180*pi 55/180*pi 3 1;
        7 50 50 70   -3/180*pi 5/180*pi 75/180*pi 3 1;
        8 55 55 20   30/180*pi 5/180*pi 55/180*pi 3 1];

    
%% Algorithm
disp('Starting backproject...');
tic

%imgs = [1 23 15 100  -3/180*pi 5/180*pi 95/180*pi 1 2];
img_pts =[1 0 0 1 1 1];
img_pts = backproject(img_pts, imgs, obj_pts, cams);
cams(:,2) = 1;

img_pts_ud = backproject(img_pts, imgs, obj_pts, cams);
cams(:,2) = 2;

toc
disp('End backproject.')

%img_pts(:,2:3) = img_pts(:,2:3) + normrnd(0, 0.0001, size(img_pts, 1), 2);
cam10 = [0.01 0 0 1e-8 1e-12 0 1e-8 1e-8];
cam20 = [0.01 0 0 1e-8 1e-12 0 1e-8 1e-8];
cams0 = [2 2 cam10 1; 3 2 cam20 1];
%cams0 = [1 1 1 1 1 2; 2 1 1 1 1 2];
%cams0 = cams;

imgs0 = [1 0 0 100   0 0 0 2 2;
         2 0 0 50    0 0 0 2 2;
         3 0 65 50   0 0 0 2 2;
         4 65 65 50  0 0 0 2 2;
         5 0 0 100   0 0 0 3 2;
         6 0 0 50    0 0 0 3 2;
         7 0 65 50   0 0 0 3 2;
         8 65 65 50  0 0 0 3 2];
    
%imgs0 = imgs;

disp('Starting bundle...');
tic

% 1. First solve without distortions
cams0(:,2) = 1;
[sol, ~] = ba_algo(img_pts, imgs0, obj_pts, cams0);
cams0(:, 3:5) = sol.cams(:, 3:5);

% 2. Solve with distortions
cams0(:,2) = 2;
[sol, stoch] = ba_algo(img_pts, sol.imgs, obj_pts, cams0);

img_pts_ud2 = backproject(img_pts, sol.imgs, obj_pts, sol.cams);
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
diagMxx = diag(stoch.Mxx);
coorMxx = sqrt(diagMxx(1:3));
rotMxx = sqrt(diagMxx(4:6));
apostMll = sqrt(diag(stoch.Mll));
fprintf('Coordinates apost. error: %.3f m\n', norm(coorMxx) );
fprintf('  Rotations apost. error: %.3f deg\n', norm(rotMxx)/pi*180 );
fprintf(' Img. coor. apost. error: %.3f pixel\n', norm(apostMll));

%
figure(3);
plot_problem(1, sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, 'g');

% Sparsity pattern
% figure(2);
show_sparisty_pattern(stoch.J, 10);

