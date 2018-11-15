clear variables; clc; close all;

KNOWN = 1;
UNKNOWN = 2;

%% Settings

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi, kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

scale = 4.87e-6;

% Camera without distortion
cam1 = [0.05 -0.1*1e-4 1e-6];
cams = [2 1 cam1 KNOWN];

% generate tie points as grid
objx = -100 : 25 : 100;
objy = -100 : 25 : 100;
[objx, objy] = meshgrid(objx, objy);
objx = objx(:);
objy = objy(:);
n_tie_pts = length(objx);

% generate control points as grid
obj_ct_x = -100 : 50: 100;
obj_ct_y = -100 : 50 : 100;
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
to_remove = find(and(abs(img_pts(:,2)/scale) > 960, abs(img_pts(:,3)/scale) > 540));
%to_remove = find(and(abs(img_pts(:,2)/scale) > 2000, abs(img_pts(:,3)/scale) > 2000));
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

%% Add error to image points
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

%% Create graph
n = size(imgs0, 1);
A = zeros(n, n);
for i  = 1 : n
    pts_i = img_pts(img_pts(:, 4) == imgs0(i, 1), :);

    for j  = i+1 : n
        pts_j = img_pts(img_pts(:, 4) == imgs0(j, 1), :);
        [~, idx_i, idx_j] = intersect(pts_i(:,5), pts_j(:,5));
        %if length(idx_i) > 0
            %fprintf('%i %i n = %i\n', i, j, length(idx_i));
            A(i, j) = length(idx_i);
            A(j, i) = length(idx_i);
        %end

    end
end
G = graph(A);
L = laplacian(G);
%plot(G,'Layout','force')

%% Clustering: Spatial clustering
n_cluster = 2^2;

minx = floor(min(imgs0(:, 2)));
maxx = ceil(max(imgs0(:, 2)));
miny = floor(min(imgs0(:, 3)));
maxy = ceil(max(imgs0(:, 3)));
dx = (maxx - minx) / sqrt(n_cluster);
dy = (maxy - miny) / sqrt(n_cluster);

clusters = {}; 
k = 0;
for i = minx : dx : maxx - dx
    for j = miny : dy : maxy - dy
        k = k + 1;
        idx = find(and( and(i <= imgs0(:, 2), imgs0(:, 2) < i + dx), and(j <= imgs0(:, 3), imgs0(:, 3) < j + dy)));
        clusters{k}.idx = idx;
        fprintf('Cluster #%2i n = %3i [%6.1f, %6.1f; %6.1f %6.1f]\n', k, length(idx), i, i+dx, j, j+dy);
    end
end

probs0.img_pts = img_pts;
probs0.imgs = imgs0;
probs0.obj_pts = obj_pts0;
probs0.cams = cams0;

for k = 1 : length(clusters)
    idx = clusters{k}.idx;
    clusters{k}.imgs    = imgs0(idx, :);
    clusters{k}.cams    = probs0.cams( ismember( probs0.cams(:, 1), clusters{k}.imgs(:, 8) ), :);
    
    ids = clusters{k}.imgs(:, 1);
    clusters{k}.img_pts = [];
    clusters{k}.obj_pts = [];    
    
    for i = 1 : length(ids)
        iidx = find(probs0.img_pts(:, 4) == ids(i));
        clusters{k}.img_pts = [clusters{k}.img_pts;   probs0.img_pts(iidx, :)];
        obj_idx = find(ismember( probs0.obj_pts(:, 1), probs0.img_pts(iidx, 5) ));
        clusters{k}.obj_pts = [clusters{k}.obj_pts; probs0.obj_pts(obj_idx, :)];        
    end    
    
    n_ctrl = 0;
    if ~isempty(clusters{k}.obj_pts)
        clusters{k}.obj_pts = unique(clusters{k}.obj_pts, 'row');
        n_ctrl = length(find((clusters{k}.obj_pts(:, end) == KNOWN)));        
    end
    
    
    fprintf('Cluster #%i n = %4i, n_imgs: %4i, n_pts: %6i, n_obj_pts: %5i, n_cams: %2i n_ctrl: %3i\n', k, length(idx), size(clusters{k}.imgs, 1), size(clusters{k}.img_pts, 1), size(clusters{k}.obj_pts, 1), size(clusters{k}.cams, 1), n_ctrl );
end
n_ctrl = length(find((probs0.obj_pts(:, end) == KNOWN)));
fprintf('Total                n_imgs: %4i, n_pts: %6i, n_obj_pts: %5i, n_cams: %2i n_ctrl: %3i\n', size(probs0.imgs, 1), size(probs0.img_pts, 1), size(probs0.obj_pts, 1), size(probs0.cams, 1), n_ctrl);

%return

%% Images 
obj_pts_per_imgs = zeros(size(probs0.imgs, 1), 2);
for k = 1 : size(probs0.imgs, 1)
    idx = probs0.imgs(k, 1);
    uobj = unique(probs0.img_pts( probs0.img_pts(:, 4) == idx, 5 ));
    obj_pts_per_imgs(k, 1) = length(uobj)*3;
    obj_pts_per_imgs(k, 2) = length(uobj)*3-9;
end
idx = find(obj_pts_per_imgs(:, 2)<0);
r = sum(abs(obj_pts_per_imgs(idx, 2)));

%%
% run_time = zeros(length(clusters), 1);
% for k = 1 : length(clusters)
%     cluster = clusters{k};
%     if ~isempty(cluster.idx)
%         tic
%         
%         [sol, stoch] = ba_algo(cluster.img_pts, cluster.imgs, cluster.obj_pts, cluster.cams);
%         J = get_jacobian(cluster.img_pts, cluster.imgs, cluster.obj_pts, cluster.cams);        
%         N = J'*J;
%         
%         run_time(k) = toc;
%         return
%     end
% end
% return

tic;
%[sol, stoch] = ba_algo(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams);
J = get_jacobian(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams); 

N = J'*J;
[~, N2] = balance(N);
n_zeig = sum(eig(N2) <= eps)

r2 = size(N2, 2) - rank(N2)
%r3 = size(N, 2) - rank(N)
r
n_deg = full(diag(L));
n_iso = sum(n_deg == 0);
%n_deg_m = n_deg - 2;
n_deg_m = sum(n_deg == 1);
n_deg_m2 = sum(n_deg == 2);
%sum(abs(n_deg_m)) + n_iso
n_iso*9 + n_deg_m*6 + n_deg_m2*3



return;

L = chol(N);
Linv = inv(L);
Ninv = Linv'*Linv;
detL = det(L)^2;

run_time_ref = toc;
return;
%% Visualization
figure(1); hold on;
opts1 = plotset;
opts1.color_cam = 'g';
%plot_problem(1, img_pts, imgs0, obj_pts0, cams0, opts1);

opts2 = plotset;
opts2.color_cam = 'r';
plot_problem(1, img_pts, imgs, obj_pts0, cams, opts2);

opts3 = plotset;
opts3.color_cam = 'b';
opts3.line_width_cam = 2;
plot_problem(1, sol.img_pts, sol.imgs, sol.obj_pts, sol.cams, opts3);





