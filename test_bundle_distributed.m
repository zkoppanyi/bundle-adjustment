clear variables; clc; close all;

KNOWN = 1;
UNKNOWN = 2;

%% Settings

% Cameras       : ID, cam_type, f, cx, cy, type
% Object points : ID, X, Y, Z, type
% Images        : ID, X, Y, Z, omega, phi, kappa, cam_id, type
% Image points  : ID, x, y, img_id, obj_id, type

scale = 4.87e-6;                                                            % Pixel sensor's size in m
overlap_ratio = 0.6;                                                        % Overlap
f = 0.05;                                                                   % Focal length
d = 6;                                                                     % In- and cross-track image distances
h = 80;                                                                     % Flight height
img_size = f * (d*overlap_ratio) / h / scale * 4;                           % Image size
pix_error = 0;

% Camera without distortion
cam1 = [f -0.1*1e-4 1e-6];
cams = [2 1 cam1 KNOWN];

% generate tie points as grid
objx = -250 : 3 : 250;
objy = -250 : 3 : 250;
[objx, objy] = meshgrid(objx, objy);
objx = objx(:) + (rand(length(objx(:)), 1)-0.5)*10;
objy = objy(:) + (rand(length(objy(:)), 1)-0.5)*10;
n_tie_pts = length(objx);

% generate control points as grid
obj_ct_x = -50 : 10 : -40;
obj_ct_y = -50 : 10 : -40;
%obj_ct_x = -200 : 50 : 200;
%obj_ct_y = -200 : 50 : 200;
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
img_pts = [1 0 0 1 1 1];
img_pts = backproject(img_pts, imgs, obj_pts, cams);
disp('End backproject.')

obj_pts0 = obj_pts;
tie_pts_idx  = 1 : n_tie_pts;
%obj_pts0(tie_pts_idx, 2:5) = repmat([0 0 0 UNKNOWN], length(tie_pts_idx), 1);
obj_pts0(tie_pts_idx, end) = UNKNOWN;
obj_pts0(tie_pts_idx, 2:4) = obj_pts(tie_pts_idx, 2:4);
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
% opts.is_show_rays = 0;
% opts.cam_scale = 1.8;
% plot_problem(2, img_pts, imgs, obj_pts0, cams, opts);
% view(-40,40)
% return

%% Add error to image points
img_pts_gt = img_pts;
img_pts(:,2:3) = img_pts(:,2:3) + normrnd(0, pix_error*scale, size(img_pts, 1), 2);

fprintf("Number of images: %i\n", size(imgs, 1));
fprintf("Number of image points: %i\n", size(img_pts, 1));
fprintf("Number of tie points: %i\n", length(tie_pts_idx));
fprintf("Number of control points: %i\n", length(control_pts_idx));

%cam10 = [0.053 -0.1*1e-4 1e-6];
%cams0 = [2 1 cam10 UNKNOWN];
cams0 = cams;
cams0(:, end) = UNKNOWN;
   
imgs0 = imgs;
imgs0(:, 2:4) = imgs0(:, 2:4) + (rand(size(imgs0,1), 3)-0.5)*2*5;
imgs0(:, 5:7) = imgs0(:, 5:7); % + (rand(size(imgs0,1), 3)-0.5)*0*5/180*pi;
imgs0(:, end) = UNKNOWN;

obj_pts0(:, 2:4) = obj_pts0(:, 2:4) + (rand(size(obj_pts0,1), 3)-0.5)*2*0;


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
n_cluster = 4^2; 

minx = floor(min(obj_pts0(:, 2)));
maxx = ceil(max(obj_pts0(:, 2)))+1;
miny = floor(min(obj_pts0(:, 3)));
maxy = ceil(max(obj_pts0(:, 3)))+1;
dx = (maxx - minx) / sqrt(n_cluster);
dy = (maxy - miny) / sqrt(n_cluster);
%dx = (maxx - minx) / 3;
%dy = (maxy - miny) / 5;

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
        clusters{k}.n_ctrl = n_ctrl;
        fprintf('Cluster #%2i n = %3i [%6.1f, %6.1f; %6.1f %6.1f] n_ctrl= %5i/%5i \n', k, length(pdix), i, i+dx, j, j+dy, n_ctrl, n_ctrl);
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

% % Check sizes
% if ( size(img_pts2, 1) ~= size(probs0.img_pts, 1) ) && ( size(img_pts2, 2) ~= size(probs0.img_pts, 2) ), error('Dimensions mismatch!'); end
% if ( size(obj_pts2, 1) ~= size(probs0.obj_pts, 1) ) && ( size(obj_pts2, 2) ~= size(probs0.obj_pts, 2) ), error('Dimensions mismatch!'); end
% if ( size(imgs2   , 1) ~= size(probs0.imgs   , 1) ) && ( size(imgs2   , 2) ~= size(probs0.imgs   , 2) ), error('Dimensions mismatch!'); end
% if ( size(cams2   , 1) ~= size(probs0.cams   , 1) ) && ( size(cams2   , 2) ~= size(probs0.cams   , 2) ), error('Dimensions mismatch!'); end
%     
% % Check indeces, if runs without problem, indexing is ok!
% J = ba_problem(img_pts2, imgs2, obj_pts2, cams2); 
% 
% % Check rank!
% N = J'*J;
% [~, N2] = balance(full(N));
% n_zeig = sum(eig(N2) <= eps); % number of zero eigenvalues
% 
% if n_zeig ~= 0
%     disp('Problem is singular!');
% end

%% 
disp('Start computing...');

probs0.img_pts  = img_pts2;
probs0.imgs     = imgs2;
probs0.obj_pts  = obj_pts2;
probs0.cams     = cams2;

%tic
%[sol, stoch]    = ba_algo(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams);
%sol.cams(2)     = 1;
%[~, ~, idxs]    = ba_problem(sol.img_pts, sol.imgs, sol.obj_pts, sol.cams); 
%x_full = idxs.x0;
%t_full = toc;

[J, r, idxs]    = ba_problem(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams); 
x0 = idxs.x0;
fprintf('norm(r)  = %.6f\n', norm(r/scale)/length(r));

x_full = x0;
r_full = r;
J_full = J;

t_distrib = 0;
t_full = 0;

for iter = 1 : 5
     
    %% Separate Jacobian 
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
        idx      = ismember(idxs.idx_imgs(:, 1), cluster.idx_imgs(:, 1));   
        idx_imgs = idxs.idx_imgs(idx, :);
        for i = 1 : size(idx_imgs, 1)
            iidx    = idx_imgs(i, 2):idx_imgs(i, 3);
            if sum(sum( J(idx_ncluster, iidx) )) ~= 0            
                cluster.coupled_imgs = [cluster.coupled_imgs, iidx];            
                %fprintf('Image %3i: Coupled\n', idxs.idx_imgs(idx(i), 1));
            else
                cluster.uncoupled_imgs = [cluster.uncoupled_imgs, iidx];                        
            end        
        end   

        % Object points
        idx         = ismember(idxs.idx_obj_pts(:, 1), cluster.idx_obj(:, 1));    
        idx_obj_pts = idxs.idx_obj_pts(idx, :);    
        for i = 1 : size(idx_obj_pts, 1)
            if idx_obj_pts(i, 2) < 1, continue; end % ok, this is a control point
            oidx = idx_obj_pts(i, 2):idx_obj_pts(i, 3);
            if sum(sum(J(idx_ncluster, oidx), 1)) ~= 0
                cluster.coupled_objs = [cluster.coupled_objs, oidx];            
                %fprintf('Object %3i: Coupled\n', idx_obj_pts(i, 1));
            else
                cluster.uncoupled_objs = [cluster.uncoupled_objs, oidx]; 
            end
        end   

        % Cameras
        idx      = find(ismember(idxs.idx_cams(:, 1), cluster.idx_cams(:, 1)));  
        idx_cams = idxs.idx_cams(idx, :);
        for i = 1 : size(idx_cams, 1)
            cidxx = idx_cams(i, 2):idx_cams(i, 3);
            if sum(sum(J(idx_ncluster, cidxx), 1)) ~= 0
                cluster.coupled_cams = [cluster.coupled_cams, cidxx];            
                %fprintf('Cam %3i: Coupled\n', idxs.idx_cams(idx(i), 1));
            else
                cluster.uncoupled_cams = [cluster.uncoupled_cams, cidxx];
            end
        end         

        coupled     = [cluster.coupled_imgs, cluster.coupled_cams, cluster.coupled_objs];
        all_coupled = [all_coupled, coupled];

        cluster.uncoupled   = [cluster.uncoupled_imgs, cluster.uncoupled_cams, cluster.uncoupled_objs];
        cluster.J   = J(idx_cluster, cluster.uncoupled);

        cluster.r   = r(idx_cluster);

        clusters{k} = cluster;

%         fprintf('\nCluster #%i\n', k);
%         fprintf('Image     : U= %5i , C= %5i\n', length(cluster.uncoupled_imgs), length(cluster.coupled_imgs))
%         fprintf('Object pt : U= %5i , C= %5i\n', length(cluster.uncoupled_objs), length(cluster.coupled_objs))
%         fprintf('Camera    : U= %5i , C= %5i\n', length(cluster.uncoupled_cams), length(cluster.coupled_cams))
          %fprintf('Cluster #%i U= %5i C= %5i n_ctrl= %i\n', k, length(cluster.uncoupled), length(coupled), cluster.n_ctrl);
          %fprintf('\\#%i & %5i & %5i &  %i \\\\ \\hline \n', k, length(cluster.uncoupled), length(coupled), cluster.n_ctrl);

    end
    all_coupled = unique(all_coupled);

    for k = 1 : length(clusters)
       cluster = clusters{k};
       idx_cluster    = cluster.idx_start : cluster.idx_end;
       clusters{k}.D  = J(idx_cluster, all_coupled);  
       clusters{k}.all_coupled = all_coupled;
    end
    
    %% Show spearated Jacobian
%     sJ = sparse(size(J, 1), size(J, 2));
%     sr = zeros(size(J, 1), 1);
%     col = 0;
%     n = size(J, 1);
%     for k = 1 : length(clusters)
%        cluster = clusters{k};
%        col_n = length(cluster.uncoupled_imgs) + length(cluster.uncoupled_cams) + length(cluster.uncoupled_objs);
%        idx_cluster    = cluster.idx_start : cluster.idx_end;
%        sJ(idx_cluster, (col+1):(col+col_n)) = cluster.J;
%        %sJ(cluster.idx_start : cluster.idx_end, cluster.uncoupled) = cluster.J;
%        sr(idx_cluster) = r(idx_cluster);
%        col = col + col_n;
%     end
%     
%     idx_start_coupled = col;
%     
%     for k = 1 : length(clusters)
%        cluster = clusters{k};
%        col_n = length(cluster.all_coupled);
%        sJ(cluster.idx_start : cluster.idx_end, (col+1):(col+col_n)) = cluster.D;
%     end
%     
%     if (idx_start_coupled+length(all_coupled)-size(J, 2)) ~= 0
%         error('Problem with indexing and sparsity discovery!')
%     end
%     
%     if col+col_n-size(J, 2) ~= 0
%         error('Problem with indexing and sparsity discovery!')
%     end
%    
%     %show_sparisty_pattern(sJ);
%     spy(sJ); hold on;
%     col = 0;
%     for k = 1 : length(clusters)
%        cluster = clusters{k};
%        col_n = length(cluster.uncoupled_imgs) + length(cluster.uncoupled_cams) + length(cluster.uncoupled_objs);   
%        
%        viz_unc = [0, cluster.idx_start; 0, cluster.idx_end; size(J, 2), cluster.idx_end; size(J, 2), cluster.idx_start; 0, cluster.idx_start];
%        plot(viz_unc(:, 1), viz_unc(:, 2), 'g.-', 'LineWidth', 3);
%     
%        viz_unc = [idx_start_coupled, cluster.idx_start; idx_start_coupled, cluster.idx_end; size(J, 2), cluster.idx_end; size(J, 2), cluster.idx_start; idx_start_coupled, cluster.idx_start];
%        plot(viz_unc(:, 1), viz_unc(:, 2), 'g.-', 'LineWidth', 3);
%        
%        viz_unc = [col, cluster.idx_start; col, cluster.idx_end; col+col_n, cluster.idx_end; col+col_n, cluster.idx_start; col, cluster.idx_start];
%        plot(viz_unc(:, 1), viz_unc(:, 2), 'r.-', 'LineWidth', 3);
%        
%        col = col + col_n;
%     end
%     
%     col = 0;
%     for k = 1 : length(clusters)
%        cluster = clusters{k};
%        col_n = length(cluster.uncoupled_imgs) + length(cluster.uncoupled_cams) + length(cluster.uncoupled_objs);                 
%        viz_unc = [col, cluster.idx_start; col, cluster.idx_end; col+col_n, cluster.idx_end; col+col_n, cluster.idx_start; col, cluster.idx_start];
%        plot(viz_unc(:, 1), viz_unc(:, 2), 'r.-', 'LineWidth', 3);       
%        col = col + col_n;
%     end

    %% Solution for separated system
%     sdx = linsolve(full(sJ'*sJ), full(sJ'*sr));
%     
%     col = 0;
%     for k = 1 : length(clusters)
%        cluster = clusters{k};
%        col_n = length(cluster.uncoupled_imgs) + length(cluster.uncoupled_cams) + length(cluster.uncoupled_objs);
%        clusters{k}.dx = sdx((col+1):(col+col_n));
%        col = col + col_n;
%     end
%         
%     col_n = length(cluster.all_coupled);
%     z = sdx((col+1):(col+col_n));
    
    %% Show problem
%     opts = plotset;
%     opts.is_show_rays = 0;
%     opts.cam_scale = 1.8;
%     cols = 'rgcmygrmgmcyrgrgcmygrmgmcyrg';
%     for k = 1 : length(clusters)
%         opts.color_cam = cols(k);     
%         ids = idxs.idx_imgs(ismember( idxs.idx_imgs(:, 2), clusters{k}.uncoupled_imgs), :);
%         cimgs = imgs( ismember(probs0.imgs(:, 1), ids(:, 1)), :);
%         plot_problem(1, [], cimgs, probs0.obj_pts, probs0.cams, opts);
%         %plot_problem(1, clusters{k}, opts);
%     end
%     opts.color_cam = 'k';
%     for k = 1 : length(clusters)
%         ids = idxs.idx_imgs(ismember( idxs.idx_imgs(:, 2), clusters{k}.coupled_imgs), :);
%         cimgs = imgs( ismember(probs0.imgs(:, 1), ids(:, 1)), :);
%         plot_problem(1, [], cimgs, probs0.obj_pts, probs0.cams, opts);
%     end
%     view(-40,40)
%     return
    
    %% Solve
    for k = 1 : length(clusters)
        clusters{k}.I = eye(size(clusters{k}.J, 2));
    end
    
    n = size(clusters{1}.D, 2);    
    for k = 1 : length(clusters)
       cluster = clusters{k};
       lJ = cluster.J;
       N = lJ'*lJ;
       D = cluster.D;
       r = cluster.r;

       % Calulcate projection: Choelsky + inverse
       
       %j   = colperm(N);
       %L   = chol(N(j, j));
       %L(j, j) = L;       
       
%        fprintf('Cholesky factorization at #%i ... ', k)
%        tic
%        L    = chol(N);      
%        invL = L \ clusters{k}.I;
%        Ninv = invL*invL';
%        P = D'*lJ*Ninv*lJ';              
%        t_distrib = t_distrib + toc;
%        fprintf('done. \n')
       
       tic
       %Ninv = N \ eye(size(N, 1));
       Ninv = inv(N);
       P = D'*lJ*Ninv*lJ';
       t_distrib = t_distrib + toc;
       
       % check result: norm(sum(Ninv - inv(N)))
       %return

       % Calulcate projection: QR + inverse
%        %tic
%        [C, R2] = qr(N, eye(size(N, 1)));
%        Ninv = R2\C;
%        %Ninv = inv(R2)*Q';
%        P = D'*lJ*Ninv*lJ';
%        %toc
       
       
       % Calulcate projection: QR with pseudo inverse.
%        %tic      
         %R2 = qr(lJ);
%        [Q, R2] = qr(lJ);
%        %toc       
%        R2 = L;
%        Q = lJ*invL;
%        %E  = eye(rank(full(R2)));
%        E  = eye(rank(full(R2)));
%        n = size(E, 1);
%        I = zeros(size(R2, 1));
%        I(1:n, 1:n) = E;
%        P = D' * Q * I * Q';    
%        %toc

       %disp('Adding...')
       tic
       clusters{k}.R = D'*D - P*D; 
       clusters{k}.l = D'*r - P*r;       
       t_distrib = t_distrib + toc;
       
       clusters{k}.Ninv = Ninv;
    end
    
    R = sparse(n, n);
    l = zeros(n, 1);
    for k = 1 : length(clusters)
        R = R + clusters{k}.R;
        l = l + clusters{k}.l;
    end
    
    disp('Calculate coupled system...')
    tic
    z = R \ l;
    %z = pinv(R)*l;
    t_distrib = t_distrib + toc;

    disp('Calculate uncoupled systems...')
    for k = 1 : length(clusters)
        lJ = clusters{k}.J;
        D = clusters{k}.D;
        r = clusters{k}.r;
        Ninv = clusters{k}.Ninv;

        tic
        clusters{k}.dx = Ninv*(lJ'*r - lJ'*D*z);        
        %clusters{k}.dx = (lJ'*lJ) \ (lJ'*r - lJ'*D*z);        
        t_distrib = t_distrib + toc;
    end

    
    
    %% Solve full problem with Gauss-Newton
    disp('Calculate full solution');
    tic;
    A = J_full'*J_full;
    b = J_full'*r_full;
    %[c, R] = qr(A, b);
    %dx = R\(R'\(A'*b))
    dx = mldivide(A, b);
    
    x_full = x_full - dx;
    [J_full, r_full]    = ba_problem(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams, x_full); 
    t_full = t_full + toc;
    
    %% Update problem
    dx = zeros(length(x0), 1);
    for k = 1 : length(clusters)
        cluster = clusters{k};    
        dx(cluster.uncoupled) = cluster.dx;
    end
    dx(all_coupled) = z;
    
    x0 = x0 - dx;
    [J, r, idxs]    = ba_problem(probs0.img_pts, probs0.imgs, probs0.obj_pts, probs0.cams, x0); 
    

    
    
    %% Report
    fprintf('Iteration #%i\n', iter);
    fprintf('norm(r)          = %10.4e    %% Norm of residuals\n', norm(r/scale)/length(r));
    fprintf('norm(dx)         = %10.4e    %% Norm of updates\n', norm(dx)); 
    fprintf('norm(x_c - x_d)  = %10.4e    %% Difference between centralized and distributed solution \n', norm(x_full - x0));
    fprintf('Time centralized : %10.3f s\n', t_full)
    fprintf('Time seq. distr. : %10.3f s  %% Time of sequential processing\n', t_distrib)
    fprintf('Time per cluster : %10.3f s  %% Time of parallel processing\n\n', t_distrib/n_cluster)
end



 






