%%
% Goal: An example and testing bundle adjusmtent for an aerial imaging
% scenario
%

%%
clear variables; clc; close all;

KNOWN = 1;
UNKNOWN = 2;

%% Settings

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

scale = 4.87e-6;

% Camera without distortion
cam1 = [0.05 -0.1*1e-4 1e-6];
cams = [2 1 cam1 KNOWN];

% generate tie points as grid
objx = -100 : 10 : 100;
objy = -100 : 10 : 100;
[objx, objy] = meshgrid(objx, objy);
objx = objx(:);
objy = objy(:);
n_tie_pts = length(objx);

% generate control points as grid
obj_ct_x = -100 : 80: 100;
obj_ct_y = -100 : 80 : 100;
[obj_ct_x, obj_ct_y] = meshgrid(obj_ct_x, obj_ct_y);
obj_ct_x = obj_ct_x(:);
obj_ct_y = obj_ct_y(:);
objx = [objx; obj_ct_x];
objy = [objy; obj_ct_y];
n_ct_pts = length(obj_ct_x);

n_test_pt = length(objx);
obj_pts = [(1:n_test_pt)', objx, objy, (rand(n_test_pt, 1)-0.5)*20, ones(n_test_pt, 1)*KNOWN];
   
imgs = [];
pt_id = 0;
for i = -50 : 10 : 50
    pt_id = pt_id + 1;
    x = i;
    y = i;
    imgs = [imgs; pt_id x y 80  -3/180*pi 5/180*pi 0/180*pi 2 UNKNOWN];
end

for i = -50 : 10 : 50
    pt_id = pt_id + 1;
    x = i;
    y = i+30;
    imgs = [imgs; pt_id x y 80  -3/180*pi 5/180*pi 0/180*pi 2 UNKNOWN];
end
    
%% Initialization with backprojection
disp('Starting backproject...');
tic
img_pts = [1 0 0 1 1 1];
img_pts = backproject(img_pts, imgs, obj_pts, cams);

toc
disp('End backproject.')


obj_pts0 = obj_pts;
tie_pts_idx  = 1 : n_tie_pts;
%obj_pts0(tie_pts_idx, 2:5) = repmat([0 0 0 UNKNOWN], length(tie_pts_idx), 1);
obj_pts0(tie_pts_idx, end) = UNKNOWN;
obj_pts0(tie_pts_idx, 2:4) = obj_pts(tie_pts_idx, 2:4)+3;
control_pts_idx = find(obj_pts0(:, end) == KNOWN);

%%  Remove points that are not on the image 
to_remove = find(and(abs(img_pts(:,2)/scale) > 2000, abs(img_pts(:,3)/scale) > 2000));
img_pts(to_remove, :) = [];

% remove points not seen by two images
rm_list = [];
for k = 1 : length(obj_pts0)
    idx = find(obj_pts0(k,1) == img_pts(:,5));
    if length(idx) < 2
        rm_list = [rm_list; idx];
    end
end
img_pts(rm_list,:) = [];

% remove object points that are not seen by any image
rm_list = [];
for k = 1 : length(obj_pts0)
    if isempty(find(obj_pts0(k,1) == img_pts(:,5)))
        rm_list = [rm_list; k];
    end
end
obj_pts0(rm_list, :) = [];
%%

% add error to image points
pix_error = 20;
img_pts_gt = img_pts;
img_pts(:,2:3) = img_pts(:,2:3) + normrnd(0, pix_error*scale, size(img_pts, 1), 2);

fprintf("Number of tie points: %i\n", length(tie_pts_idx));
fprintf("Number of control points: %i\n", length(control_pts_idx));

cam10 = [0.053 -0.1*1e-4 1e-6];
cams0 = [2 1 cam10 UNKNOWN];
   
imgs0 = imgs;
imgs0(:, 2:4) = imgs0(:, 2:4) + (rand(size(imgs0,1), 3)-0.5)*2*5;
imgs0(:, 5:7) = imgs0(:, 5:7) + (rand(size(imgs0,1), 3)-0.5)*2*5/180*pi;
imgs0(:, end) = UNKNOWN;

disp('Starting bundle...');
tic

[sol, stoch] = ba_algo(img_pts, imgs0, obj_pts0, cams0);

x0 = ba_problem(img_pts, imgs0, obj_pts0, cams0, 2);
r = ba_problem(img_pts, imgs0, obj_pts0, cams0, x0, 1);
J = ba_problem(img_pts, imgs0, obj_pts0, cams0, x0, 2);
stoch.J = J;

opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true, 'Display', 'iter');
%opts.Algorithm = 'levenberg-marquardt';
x2 = lsqnonlin(@(x) ba_problem(img_pts, imgs0, obj_pts0, cams0, x, 1) ,x0, [], [], opts);

sol.cams(2) = 1;
x = ba_problem(sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, 2);


%%

N = stoch.J'*stoch.J;
if size(N,1) ~= rank(N)
    disp("Wrong J structure!")
else
    disp("J structure is fine!")    
end

toc
disp('End bundle.')

%imgs
%imgs0
%sol.imgs
dp = sol.imgs(:,2:4) - imgs(:,2:4);
fprintf("Image coordinate error from ground truth: %.3f\n", mean(abs(dp(:))) );
%dp = sol.img_pts(:,2:3) - img_pts_gt(:,2:3);
%fprintf("Simulated pixel error: %.3f\n", pix_error);
%fprintf("Predicted pixel error: %.3f\n", std(dp(:))/scale);


%return;
%% Visualization

% Stochastics results
% diagMxx = diag(stoch.Mxx);
% coorMxx = sqrt(diagMxx(1:3));
% rotMxx = sqrt(diagMxx(4:6));
% apostMll = sqrt(diag(stoch.Mll));
% fprintf('Coordinates apost. error: %.3f m\n', norm(coorMxx) );
% fprintf('  Rotations apost. error: %.3f deg\n', norm(rotMxx)/pi*180 );
% fprintf(' Img. coor. apost. error: %.3f pixel\n', norm(apostMll));

%
figure(1); hold on;
opts1 = plotset;
opts1.color_cam = 'g';
plot_problem(1, img_pts, imgs0, obj_pts0, cams0, opts1);

opts2 = plotset;
opts2.color_cam = 'r';
plot_problem(1, img_pts, imgs, obj_pts0, cams, opts2);

opts3 = plotset;
opts3.color_cam = 'b';
opts3.line_width_cam = 2;
plot_problem(1, sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, opts3);

% Sparsity pattern
show_sparisty_pattern(stoch.J, 2);
show_sparisty_pattern(stoch.J'*stoch.J, 3);

N = stoch.J'*stoch.J;
if size(N,1) ~= rank(N)
    disp("Wrong J structure!")
end

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


