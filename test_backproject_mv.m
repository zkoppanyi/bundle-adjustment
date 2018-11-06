%%
% Goal: Backprojection for multiview cameras
%
%%

clear variables; clc; close all;

scale = 4.87e-6;

KNOWN = 1;
UNKNOWN = 2;

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

% Cameras
% ID, cam_type, f, cx, cy, type
%cam1 = [1 -0.1*1e-4 1e-6];
%cam2 = [1  0.1*1e-4 1e-6];
cam1 = [1  0 0];
cam2 = [1  0 0];
cams = [1 1 cam1 KNOWN; 
        2 1 cam2 KNOWN];


% Object points
% ID, X, Y, Z, type
n_test_pt = 15;
rand_pts = [(rand(n_test_pt, 1)-0.5)*150, (rand(n_test_pt, 1)-0.5)*150, (rand(n_test_pt, 1)-0.5)*30];
obj_pts = [(1:n_test_pt)' rand_pts ones(n_test_pt, 1)];

%Images 
% ID, X, Y, Z, omega, phi kappa, cam_id, type
imgs = [1 0 0 100   0/180*pi 0/180*pi 0/180*pi 1 1;
        2 20 50 100 0/180*pi 0/180*pi 6/180*pi 2 1];
b_vec = imgs(2, 2:4) - imgs(1, 2:4);
b_len = norm(b_vec);
d_ang = imgs(2, 5:7) - imgs(1, 5:7);

% Image points
% ID, x, y, img_id, obj_id, type
img_pts = [1 0 0 1 1 1];

img_pts = backproject(img_pts, imgs, obj_pts, cams);

plot_problem(1, img_pts, imgs, obj_pts, cams);

%% Calculate fundamental matrix for checking the result
width = 10; height = 10;
img_pts1 = img_pts(img_pts(:,4) == 1, :);
img_pts2 = img_pts(img_pts(:,4) == 2, :);
img_pts1(:,2:3) = (img_pts1(:,2:3) - cams(1,4:5)) + [width height];
img_pts2(:,2:3) = (img_pts2(:,2:3) - cams(2,4:5)) + [width height];
F = estimateFundamentalMatrix(img_pts1(:,2:3), img_pts2(:,2:3), 'Method', 'Norm8Point');

%K1 = [cams(1, 3) 0 0; 0 cams(1, 3) 0; w h 1];
%K2 = [cams(2, 3) 0 0; 0 cams(2, 3) 0; w h 1];

[U, ~, V] = svd(F,'econ');
W = [0, -1, 0; 1, 0, 0; 0, 0, 1];
R1 = U*W*V';
R2 = U*W'*V';
rotm2eul(R1,'ZYX')/pi*180
rotm2eul(R2,'ZYX')/pi*180

% Decompose the matrix E by svd
[u, s, v] = svd(F);

%
w = [0 -1 0; 1 0 0; 0 0 1];
z = [0 1 0; -1 0 0; 0 0 1];

% 
% E = SR where S = [t]_x and R is the rotation matrix.
% E can be factorized as:
%s = u * z * u';

% Two possibilities:
rot1 = u * w  * v'; rotm2eul(rot1,'ZYX') / pi*180
rot2 = u * w' * v'; rotm2eul(rot2,'ZYX') / pi*180

% % Two possibilities:
% t1 = u(:,3) ./max(abs(u(:,3)));
% t2 = -u(:,3) ./max(abs(u(:,3)));
% 
% 
% % 4 possible choices of the camera matrix P2 based on the 2 possible
% % choices of R and 2 possible signs of t.
% rot(:,:,1) = rot1; 
% t(:,:,1) = t1;
% 
% rot(:,:,2) = rot2; 
% t(:,:,2) = t2;
% 
% rot(:,:,3) = rot1; 
% t(:,:,3) = t2;
% 
% rot(:,:,4) = rot2; 
% t(:,:,4) = t1;

% get relative camera pose
K1 = cameraIntrinsics(cam1(1),cam1(2:3) + [width height], [2*width 2*height]);
K2 = cameraIntrinsics(cam2(1),cam2(2:3) + [width height], [2*width 2*height]);
[rot, t] = relativeCameraPose(F, K1, K2, img_pts1(:,2:3), img_pts2(:,2:3));
eul_ang = rotm2eul(rot,'XYZ');
eul_ang(1:2) = -eul_ang(1:2);

% Transform the baseline to the global system
R = eul2rotm(imgs(1, 5:7), 'XYZ');
tg = R'*t';

fprintf("Euler angles from F matrix [deg]: %.3f %.3f %.3f Errors: %.3f %.3f %.3f\n", eul_ang/pi*180, angdiff(eul_ang, d_ang)/pi*180);
fprintf("      Baseline from F matrix [m]: %.3f %.3f %.3f Errors: %.3f %.3f %.3f\n", tg*b_len, tg*b_len - b_vec');

if norm(angdiff(eul_ang, d_ang)) > 10e-5
    error("Angles are not correct from fundmental matrix!")
end

if norm(tg*b_len - b_vec') > 10e-5
    error("baseline is not correct from fundmental matrix!")
end
