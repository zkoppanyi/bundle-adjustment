%%
% Goal: Testing the bundle adjustment on large scale problem
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
for k = -100 : 15 : 100 % cross flight
    for i = -100 : 15 : 100
        pt_id = pt_id + 1;
        x = i;
        y = i + k;
        imgs = [imgs; pt_id x y 80  -3/180*pi 5/180*pi 0/180*pi 2 UNKNOWN];
    end
end
imgs(:, 4) = imgs(:, 4) + (rand(size(imgs, 1), 1)-0.5)*1;
    
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
%obj_pts0(tie_pts_idx, 2:4) = obj_pts(tie_pts_idx, 2:4)+3;
obj_pts0(tie_pts_idx, 2:4) = obj_pts(tie_pts_idx, 2:4) + (rand(length(tie_pts_idx), 3)-0.5)*5;
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

% Remove object points that are not seen by any image
rm_list = [];
for k = 1 : length(obj_pts0)
    if isempty(find(obj_pts0(k,1) == img_pts(:,5)))
        rm_list = [rm_list; k];
    end
end
obj_pts0(rm_list, :) = [];

%%

% add error to image points
pix_error = 100;
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

% x0 = ba_problem(img_pts, imgs0, obj_pts0, cams0, 2);
% %r = ba_problem(img_pts, imgs0, obj_pts0, cams0, x0, 1);
% %J = ba_problem(img_pts, imgs0, obj_pts0, cams0, x0, 2);
% %stoch.J = J;
% 
% opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true, 'Display', 'iter');
% %opts.Algorithm = 'levenberg-marquardt';
% x2 = lsqnonlin(@(x) ba_problem(img_pts, imgs0, obj_pts0, cams0, x, 1) ,x0, [], [], opts);
% sol.cams(2) = 1;
% x = ba_problem(sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, 2);
% norm(x-x2)

disp('End bundle.')



%% Errors

N = stoch.J'*stoch.J;
if find(sum(stoch.J) == 0)
    disp("Wrong J structure!")
end
% if size(N,1) ~= rank(N)
%     disp("Wrong J structure!")
% end

dxyz = obj_pts - sol.obj_pts;
dxyz0 = obj_pts - obj_pts0;
fprintf("Avg. error of obj_pts: %8.3f %8.3f %8.3f [%8.3f %8.3f %8.3f]\n", mean(dxyz(:,2)), mean(dxyz(:,3)), mean(dxyz(:,4)),                mean(dxyz0(:,2)), mean(dxyz0(:,3)), mean(dxyz0(:,4)));
fprintf("Abs. error of obj_pts: %8.3f %8.3f %8.3f [%8.3f %8.3f %8.3f]\n", mean(abs(dxyz(:,2))), mean(abs(dxyz(:,3))), mean(abs(dxyz(:,4))), mean(abs(dxyz0(:,2))), mean(abs(dxyz0(:,3))), mean(abs(dxyz0(:,4))));
fprintf("Std. error of obj_pts: %8.3f %8.3f %8.3f [%8.3f %8.3f %8.3f]\n", std(dxyz(:,2)), std(dxyz(:,3)), std(dxyz(:,4)),                std(dxyz0(:,2)), std(dxyz0(:,3)), std(dxyz0(:,4)));

dpix = (img_pts_gt(:,2:3) - sol.img_pts(:,2:3))/scale;
dpix0 = (img_pts_gt(:,2:3) - img_pts(:,2:3))/scale;
fprintf("Avg. error of pixels : %8.1f %8.1f [%8.1f %8.1f]\n", mean(dpix(:,1)), mean(dpix(:,2)),           mean(dpix0(:,1)), mean(dpix0(:,2))  );
fprintf("Abs. error of pixels : %8.1f %8.1f [%8.1f %8.1f]\n", mean(abs(dpix(:,1))), mean(abs(dpix(:,2))), mean(abs(dpix0(:,1))), mean(abs(dpix0(:,2)))  );
fprintf("Std. error of pixels : %8.1f %8.1f [%8.1f %8.1f]\n", std(dpix(:,1)), std(dpix(:,2)),           std(dpix0(:,1)), std(dpix0(:,2))  );

%return;
%% Visualization
figure(1); clf; hold on;
opts1 = plotset;
opts1.cam_scale = 5;
opts1.color_cam = 'g';
%plot_problem(1, img_pts, imgs0, obj_pts0, cams0, opts1);

opts2 = opts1;
opts2.color_cam = 'r';
plot_problem(1, img_pts, imgs, obj_pts0, cams, opts2);

opts3 = opts1;
opts3.color_cam = 'b';
opts3.line_width_cam = 2;
plot_problem(1, sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, opts3);

%show_sparisty_pattern(stoch.J, 2);
%show_sparisty_pattern(stoch.J'*stoch.J, 3);





