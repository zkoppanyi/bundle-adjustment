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
img_size = f * (d*overlap_ratio) / h / scale * 3;                           % Image size

% Camera without distortion
cam1 = [f -0.1*1e-4 1e-6];
cams = [2 1 cam1 KNOWN];

% generate tie points as grid
objx = -250 : 10 : 250;
objy = -250 : 10 : 250;
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

% %plot(G,'Layout','force')

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





