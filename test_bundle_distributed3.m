clear variables; clc; close all;

KNOWN = 1;
UNKNOWN = 2;

%% Settings

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi, kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

scale = 4.87e-6;                                                            % Pixel sensor's size in m
overlap_ratio = 0.8;                                                        % Overlap
f = 0.05;                                                                   % Focal length
d = 15;                                                                     % In- and cross-track image distances
h = 80;                                                                     % Flight height
img_size = f * (d*overlap_ratio) / h / scale * 4;                           % Image size

% Camera without distortion
cam1 = [f -0.1*1e-4 1e-6];
cams = [2 1 cam1 KNOWN];

% generate tie points as grid
objx = -250 : 15 : 250;
objy = -250 : 15 : 250;
[objx, objy] = meshgrid(objx, objy);
objx = objx(:); % + (rand(length(objx(:)), 1)-0.5)*1;
objy = objy(:); % + (rand(length(objy(:)), 1)-0.5)*1;
n_tie_pts = length(objx);

% generate control points as grid
obj_ct_x = -200 : 50 : 200;
obj_ct_y = -200 : 50 : 200;
[obj_ct_x, obj_ct_y] = meshgrid(obj_ct_x, obj_ct_y);
obj_ct_x = obj_ct_x(:) + (rand(length(obj_ct_x(:)), 1)-0.5)*10;
obj_ct_y = obj_ct_y(:) + (rand(length(obj_ct_y(:)), 1)-0.5)*10;
objx = [objx; obj_ct_x];
objy = [objy; obj_ct_y];
n_ct_pts = length(obj_ct_x);

n_test_pt = length(objx);
obj_pts = [(1:n_test_pt)', objx, objy, (rand(n_test_pt, 1)-0.5)*20, ones(n_test_pt, 1)*KNOWN];
   
imgs = [];
pt_id = 0;
for k = -100 : d : 100 % cross flight
    for i = -100 : d : 100
        pt_id = pt_id + 1;
        x = i       + (rand-0.5)*1;
        y = k       + (rand-0.5)*1;
        %y = i + k       + (rand-0.5)*1;
        z = 80      + (rand-0.5)*1;
        %imgs = [imgs; pt_id x y z  0/180*pi 0/180*pi 0/180*pi 2 UNKNOWN];
        imgs    = [imgs; pt_id x y z  -3/180*pi 5/180*pi 0/180*pi 2 UNKNOWN];
        %a = (rand(3,1)-0.5)*1/180*pi;
        %imgs = [imgs; pt_id x y z  a(1) a(2) a(3) 2 UNKNOWN];
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
to_keep = find( and( abs(img_pts(:,2)/scale) < img_size/2, abs(img_pts(:,3)/scale) < img_size/2) );
img_pts = img_pts(to_keep, :);

% remove points not seen by two images
rm_list_pts = [];
rm_list_obj = [];
for k = 1 : length(obj_pts0)
    idx = find(obj_pts0(k,1) == img_pts(:,5));
    if length(idx) < 2
        rm_list_pts = [rm_list_pts; idx];
        rm_list_obj = [rm_list_obj; k];
    end
end
img_pts(rm_list_pts,:) = [];
obj_pts0(rm_list_obj, :) = [];

% % Show the whole problem
% figure(2); hold on;
% opts = plotset;
% opts.is_show_rays = 1;
% plot_problem(2, img_pts, imgs, obj_pts0, cams, opts);
% return

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
probs0.img_pts = img_pts;
probs0.imgs = imgs0;
probs0.obj_pts = obj_pts0;
probs0.cams = cams0;

% G = create_img_graph(probs0);
% L = laplacian(G);
% if sum(diag(L)==0) ~= 0
%     error('Isolated image!');    
% end
%plot(G,'Layout','force')

%% Clustering: Spatial clustering
n_cluster = 2^2;

minx = floor(min(obj_pts0(:, 2)));
maxx = ceil(max(obj_pts0(:, 2)))+1;
miny = floor(min(obj_pts0(:, 3)));
maxy = ceil(max(obj_pts0(:, 3)))+1;
%dx = (maxx - minx) / sqrt(n_cluster);
%dy = (maxy - miny) / sqrt(n_cluster);
dx = (maxx - minx) / 1;
dy = (maxy - miny) / n_cluster;

clusters = {}; 
k = 0;
idx_start = 0;
for i = minx : dx : maxx - dx
    for j = miny : dy : maxy - dy
        k = k + 1;
        
        oidx = find(and( and(i <= obj_pts0(:, 2), obj_pts0(:, 2) < i + dx), and(j <= obj_pts0(:, 3), obj_pts0(:, 3) < j + dy)));        
        pdix = find(ismember(img_pts(:, 5), obj_pts0(oidx, 1)) );        
        iidx = find(ismember(probs0.imgs(:, 1), img_pts(pdix, 4)));
        cidx = find(ismember(probs0.cams(:, 1), probs0.imgs(iidx, 8) ));
        
        % Unknowns' locations
        clusters{k}.idx_obj     = [obj_pts0(oidx, 1),    oidx];                
        clusters{k}.idx_img_pts = [img_pts(pdix, 1),     pdix];        
        clusters{k}.idx_imgs    = [probs0.imgs(iidx, 1), iidx];       
        clusters{k}.idx_cams    = [probs0.cams(cidx, 1), cidx];
        
        % Measurements' locations
        clusters{k}.idx_start   = idx_start + 1;
        clusters{k}.idx_end     = idx_start + length(pdix)*2;
        idx_start               = idx_start + length(pdix)*2;        
        
        n_ctrl = length(find(obj_pts0(oidx, end) == KNOWN));        
        n_sum_ctrl = length(find(obj_pts0(:, end) == KNOWN));        
        fprintf('Cluster #%2i n = %3i [%6.1f, %6.1f; %6.1f %6.1f] n_ctrl= %5i/%5i \n', k, length(pdix), i, i+dx, j, j+dy, n_ctrl, n_sum_ctrl);
    end
end

% Measurements (rows) were reordered, so create a new problem
img_pts2 = []; obj_pts2 = []; imgs2 = []; cams2 = [];
for k = 1 : length(clusters)
        cluster    = clusters{k};
        img_pts2   = [img_pts2; probs0.img_pts(cluster.idx_img_pts(:, 2), :)];
        obj_pts2   = [obj_pts2; probs0.obj_pts(cluster.idx_obj(:, 2), :)];
        imgs2      = [imgs2;    probs0.imgs(cluster.idx_imgs(:, 2), :)];
        cams2      = [cams2;    probs0.cams];
        
        cluster.img_pts  = probs0.img_pts(cluster.idx_img_pts(:, 2), :);
        cluster.obj_pts  = probs0.obj_pts(cluster.idx_obj(:, 2), :);
        cluster.imgs     = probs0.imgs(cluster.idx_imgs(:, 2), :);
        cluster.cams     = probs0.cams;
        clusters{k}      = cluster;
end
imgs2 	 = unique(imgs2, 'row');
obj_pts2 = unique(obj_pts2, 'row');
cams2    = unique(cams2, 'row');

%% Check that the problem is same after clustering

% Check sizes
if ( size(img_pts2, 1) ~= size(probs0.img_pts, 1) ) && ( size(img_pts2, 2) ~= size(probs0.img_pts, 2) ), error('Dimensions mismatch!'); end
if ( size(obj_pts2, 1) ~= size(probs0.obj_pts, 1) ) && ( size(obj_pts2, 2) ~= size(probs0.obj_pts, 2) ), error('Dimensions mismatch!'); end
if ( size(imgs2   , 1) ~= size(probs0.imgs   , 1) ) && ( size(imgs2   , 2) ~= size(probs0.imgs   , 2) ), error('Dimensions mismatch!'); end
if ( size(cams2   , 1) ~= size(probs0.cams   , 1) ) && ( size(cams2   , 2) ~= size(probs0.cams   , 2) ), error('Dimensions mismatch!'); end
    
% Check indeces, if runs without problem, indexing is ok!
J = ba_problem(img_pts2, imgs2, obj_pts2, cams2); 

% Check rank!
N = J'*J;
[~, N2] = balance(full(N));
n_zeig = sum(eig(N2) <= eps); % number of zero eigenvalues

if n_zeig ~= 0
    disp('Problem is singular!');
end

%% Separate Jacobian 
probs0.img_pts  = img_pts2;
probs0.imgs     = imgs2;
probs0.obj_pts  = obj_pts2;
probs0.cams     = cams2;
[J, r, idxs]    = ba_problem(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams); 

n = 1 : size(J, 1);
%for k = 1 : length(clusters)
all_coupled = [];
for k = 1 : length(clusters)
    
    cluster = clusters{k};
    cluster.coupled_imgs = []; cluster.uncoupled_imgs = []; 
    cluster.coupled_objs = []; cluster.uncoupled_objs = []; 
    cluster.coupled_cams = []; cluster.uncoupled_cams = []; 
    
    idx_cluster    = cluster.idx_start : cluster.idx_end;
    idx_ncluster   = setdiff(n, idx_cluster);
    
    % Images
    idx     = find( ismember(idxs.idx_imgs(:, 1), cluster.idx_imgs(:, 1) ));   
    for i = 1 : length(idx)
        iidx    = idxs.idx_imgs(idx(i), 2):idxs.idx_imgs(idx(i), 3);
        if sum(sum( J(idx_ncluster, iidx) )) ~= 0            
            cluster.coupled_imgs = [cluster.coupled_imgs, iidx];            
            %fprintf('Image %3i: Coupled\n', idxs.idx_imgs(idx(i), 1));
        else
            cluster.uncoupled_imgs = [cluster.uncoupled_imgs, iidx];                        
        end        
    end   
    
    % Object points
    idx         = ismember(idxs.idx_obj_pts(:, 1), cluster.idx_obj(:, 1) );    
    idx_obj_pts = idxs.idx_obj_pts(idxs.idx_obj_pts(idx, 2) ~= -1, :);    % remove control points
    for i = 1 : size(idx_obj_pts, 1)
        oidx = idx_obj_pts(i, 2):idx_obj_pts(i, 3);
        if sum(sum(J(idx_ncluster, oidx), 1)) ~= 0
            cluster.coupled_objs = [cluster.coupled_objs, oidx];            
            %fprintf('Object %3i: Coupled\n', idx_obj_pts(i, 1));
        else
             cluster.uncoupled_objs = [cluster.uncoupled_objs, oidx]; 
        end
    end   
    
    % Cameras
    idx   = find(ismember(idxs.idx_cams(:, 1), cluster.idx_cams(:, 1) ));          
    for i = 1 : size(idx, 1)
        cidxx = idxs.idx_cams(idx(i), 2):idxs.idx_cams(idx(i), 3);
        if sum(sum(J(idx_ncluster, cidxx), 1)) ~= 0
            cluster.coupled_cams = [cluster.coupled_cams, cidxx];            
            %fprintf('Cam %3i: Coupled\n', idxs.idx_cams(idx(i), 1));
        else
            cluster.uncoupled_cams = [cluster.uncoupled_cams, cidxx];
        end
    end     
    
    coupled     = [cluster.coupled_imgs, cluster.coupled_cams, cluster.coupled_objs];
    all_coupled = [all_coupled, coupled];
    
    uncoupled = [cluster.uncoupled_imgs, cluster.uncoupled_cams, cluster.uncoupled_objs];
    cluster.J = J(idx_cluster, uncoupled);
    
    cluster.r = r(idx_cluster);
    
    clusters{k} = cluster;
    
    fprintf('\nCluster #%i\n', k);
    fprintf('Image     : U= %5i , C= %5i\n', length(cluster.uncoupled_imgs), length(cluster.coupled_imgs))
    fprintf('Object pt : U= %5i , C= %5i\n', length(cluster.uncoupled_objs), length(cluster.coupled_objs))
    fprintf('Camera    : U= %5i , C= %5i\n', length(cluster.uncoupled_cams), length(cluster.coupled_cams))
end
all_coupled = unique(all_coupled);

for k = 1 : length(clusters)
   cluster = clusters{k};
   idx_cluster    = cluster.idx_start : cluster.idx_end;
   clusters{k}.D = J(idx_cluster, all_coupled);  
   clusters{k}.all_coupled = all_coupled;
end

%% Show spearated Jacobian
sJ = sparse(size(J, 1), size(J, 2));
col = 0;
n = size(J, 1);
for k = 1 : length(clusters)
   cluster = clusters{k};
   col_n = length(cluster.uncoupled_imgs) + length(cluster.uncoupled_cams) + length(cluster.uncoupled_objs);
   sJ(cluster.idx_start : cluster.idx_end, (col+1):(col+col_n)) = cluster.J;
   col = col + col_n;
end

idx_start_coupled = col;

for k = 1 : length(clusters)
   cluster = clusters{k};
   col_n = length(cluster.all_coupled);
   sJ(cluster.idx_start : cluster.idx_end, (col+1):(col+col_n)) = cluster.D;
end

idx_start_coupled+length(all_coupled)-size(J, 2) % ?= 0
col+col_n-size(J, 2)  % ?= 0
return

show_sparisty_pattern(sJ);
col = 0;
for k = 1 : length(clusters)
   cluster = clusters{k};
   col_n = length(cluster.uncoupled_imgs) + length(cluster.uncoupled_cams) + length(cluster.uncoupled_objs);   
   
   viz_unc = [0, n-cluster.idx_start; 0, n-cluster.idx_end; size(J, 2), n-cluster.idx_end; size(J, 2), n-cluster.idx_start; 0, n-cluster.idx_start];
   plot(viz_unc(:, 1), viz_unc(:, 2), 'g.-', 'LineWidth', 3);

   viz_unc = [idx_start_coupled, n-cluster.idx_start; idx_start_coupled, n-cluster.idx_end; size(J, 2), n-cluster.idx_end; size(J, 2), n-cluster.idx_start; idx_start_coupled, n-cluster.idx_start];
   plot(viz_unc(:, 1), viz_unc(:, 2), 'g.-', 'LineWidth', 3);
   
   viz_unc = [col, n-cluster.idx_start; col, n-cluster.idx_end; col+col_n, n-cluster.idx_end; col+col_n, n-cluster.idx_start; col, n-cluster.idx_start];
   plot(viz_unc(:, 1), viz_unc(:, 2), 'r.-', 'LineWidth', 3);
   
   col = col + col_n;
end

%% Show problem
% opts = plotset;
% opts.is_show_rays = 0;
% cols = 'rgbcyk';
% for k = 1 : length(clusters)
%     opts.color_cam = cols(k);
%     plot_problem(1, clusters{k}, opts);
% end
% opts.color_cam = 'k';
% plot_problem(1, [], server.imgs, probs0.obj_pts, probs0.cams, opts);

%% Solve
lJ = J;
n = size(clusters{1}.D, 2);
R = zeros(n, n);
l = zeros(n, 1);
for k = 1 : length(clusters)
   cluster = clusters{k};
   J = cluster.J;
   N = J'*J;
   D = cluster.D;
   r = cluster.r;
   
   % calulcate inverse: Choelsky
   L    = chol(N);
   invL = inv(L);
   Ninv = invL*invL';
   % check result: norm(sum(Ninv - inv(N)))
   
   A = D'*J*Ninv*J';
   R = R + D'*D - A*D; 
   l = l + D'*r - A*r;
   
   %figure(k);
   %spy(R)
end
%z = linsolve(R, l);
z = pinv(R)*l;

return;


%% 
% Create clusters
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

for i = 1 : length(clusters)
    servers{i}.img_pts   = [];
    servers{i}.imgs      = [];
    servers{i}.obj_pts   = [];
    servers{i}.cams      = [];
end

% Create connections 
for i = 1 : length(clusters)
    obj_i = clusters{i}.obj_pts;

    for j = i+1 : length(clusters)
        obj_j = clusters{j}.obj_pts;
        [ids, iidx, jidx] = intersect(obj_i(:,1), obj_j(:,1));
        
        servers{i}.obj_pts = [servers{i}.obj_pts; obj_i(iidx,:)];
        servers{j}.obj_pts = [servers{j}.obj_pts; obj_j(jidx,:)];

        % Remove object points
        obj_i(iidx,:) = [];
        obj_j(jidx,:) = [];
        
        % Remove image points corresponding to this object point
        for k = 1 : length(ids)            
            idx = find(clusters{i}.img_pts(:, 5) == ids(k) );                
            servers{i}.img_pts = [servers{i}.img_pts; clusters{i}.img_pts(idx, :)];
            clusters{i}.img_pts(idx, :) = [];            
            
            idx = find(clusters{j}.img_pts(:, 5) == ids(k) );                
            servers{j}.img_pts = [servers{j}.img_pts; clusters{j}.img_pts(idx, :)];
            clusters{j}.img_pts(idx, :) = [];            
        end   
        
        clusters{j}.obj_pts = obj_j;
    end
    clusters{i}.obj_pts = obj_i;
    
end

for i = 1 : length(clusters)
    
     while 1 
         
            % remove images that do not have at least 5 points
            rm_list_imgs = [];
            rm_list_imgs_pts = [];
            for k = 1 : size(clusters{i}.imgs, 1)
                idx = find(clusters{i}.img_pts(:, 4) == clusters{i}.imgs(k, 1));
                if length(idx) < 5
                    rm_list_imgs = [rm_list_imgs; k];
                    rm_list_imgs_pts = [rm_list_imgs_pts; idx];
                end
            end
            
            % add image and images and objsect points to the server
            servers{i}.imgs     = [servers{i}.imgs; clusters{i}.imgs(rm_list_imgs, :)];
            servers{i}.img_pts  = [servers{i}.img_pts; clusters{i}.img_pts(rm_list_imgs_pts, :)];
            
            % ok now we are ready to remove images and image points
            clusters{i}.imgs(rm_list_imgs, :)  = [];
            clusters{i}.img_pts(rm_list_imgs_pts, :)  = [];

            % Remove points not seen by two images
            rm_list_pts = [];
            rm_list_obj = [];
            for k = 1 : length(clusters{i}.obj_pts)
                idx = find(clusters{i}.obj_pts(k,1) == clusters{i}.img_pts(:,5));
                if length(idx) < 2
                    rm_list_obj = [rm_list_obj; k];
                    rm_list_pts = [rm_list_pts; idx];
                end
            end
            
             % add image and object points to the server
            servers{i}.img_pts = [servers{i}.img_pts; clusters{i}.img_pts(rm_list_pts, :)];
            servers{i}.obj_pts = [servers{i}.obj_pts; clusters{i}.obj_pts(rm_list_obj, :)];
            
            % ok now we are ready to remove image and object points
            clusters{i}.img_pts(rm_list_pts, :)  = [];
            clusters{i}.obj_pts(rm_list_obj, :)  = [];
            
            if isempty(rm_list_imgs) && isempty(rm_list_imgs_pts) && isempty(rm_list_pts) && isempty(rm_list_obj), break; end
     end
     
end

for i = 1 : length(clusters)
    servers{i}.img_pts   = unique(servers{i}.img_pts  , 'row');
    
    idx = find(ismember( imgs0(:, 1), servers{i}.img_pts(:, 4) ));
    servers{i}.imgs      = [servers{i}.imgs; imgs0(idx, :)];
    servers{i}.imgs      = unique(servers{i}.imgs     , 'row');
    
    idx = find(ismember( obj_pts0(:, 1), servers{i}.img_pts(:, 5) ));
    servers{i}.obj_pts   = [servers{i}.obj_pts; obj_pts0(idx, :)];
    servers{i}.obj_pts   = unique(servers{i}.obj_pts  , 'row');

    servers{i}.img_pts   = clusters{i}.img_pts;

    %servers{i}.cams      = unique(servers{i}.cams     , 'row');
    servers{i}.cams      = clusters{i}.cams;
end


%return


%% Calculation

% Check that the problem is full rank
% J = get_jacobian(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams); 
% N = J'*J;
% n_zeig = sum(eig(N) <= eps);
% if n_zeig ~= 0, error('Problem is not full rank!'); end


run_time = zeros(length(clusters), 1);
cols = 'rgbcyk';
for k = 1 : length(clusters)
    cluster = clusters{k};
    server  = servers{k};
    
    if ~isempty(cluster.idx)
        
        figure(1); hold on;
        opts = plotset;
        opts.color_cam = cols(k);
        %opts.is_show_rays = 1;
        plot_problem(1, cluster.img_pts, cluster.imgs, cluster.obj_pts, cluster.cams, opts);
        opts.color_cam = 'k';
        plot_problem(1, server.img_pts, server.imgs, server.obj_pts, server.cams, opts);
        
        G = create_img_graph(cluster);
        L = laplacian(G);
        if sum(diag(L)<2) ~= 0
            fprintf('Isolated image or not enough tie points in cluster #%i!\n', k);    
        end

        tic
        
        %[sol, stoch] = ba_algo(cluster.img_pts, cluster.imgs, cluster.obj_pts, cluster.cams);
        J = get_jacobian(cluster.img_pts, cluster.imgs, cluster.obj_pts, cluster.cams); 
        D = get_jacobian(server.img_pts, server.imgs, server.obj_pts, server.cams);
        
        N = J'*J;
        [~, N2] = balance(full(N));
        n_zeig = sum(eig(N2) <= eps);
        
        if n_zeig > 0
            n_zeig
            error('Matrix is not positive definite or unbalanced in cluster #%i', k);
        end
        
        L    = chol(N2);
%         invL = inv(L);
%         Ninv = invL'*invL;
        
        run_time(k) = toc;
        
       
        %return
    end
end
return
% tic;
% [sol, stoch] = ba_algo(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams);
% J = get_jacobian(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams); 
% run_time_ref = toc;

return
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





