%%
% Goal: Test triangulation (TODO: the interface has been changed, so it doesn't
% work currently; needs to be fixed)
%%

clear variables; clc;

%% Stress triangulation

%cam1 = [0.5 1 2 0 0 100 3/180*pi 5/180*pi 25/180*pi];
%cam2 = [0.5 1 2 50 0 100 -3/180*pi -8/180*pi -30/180*pi];
%cam3 = [0.5 1 -5 50 50 100 -3/180*pi -8/180*pi -15/180*pi];
%cams = [cam1; cam2; cam3];

cams = [];
for i = 1 : 5,
    cam = [0.5 1 2 (rand-0.5)*100 (rand-0.5)*100 100-(rand-0.5)*10 (rand-0.5)*45/180*pi (rand-0.5)*45/180*pi (rand-0.5)*45/180*pi];
    cams = [cams; cam];
end

n_test_pt = 1;
pt = [(rand(n_test_pt, 1)-0.5)*150, (rand(n_test_pt, 1)-0.5)*150, (rand(n_test_pt, 1)-0.5)*30];
%pt = [0, 0, 50];


%plot_problem(1, cams, pts, pti);

% cameras
cam_ids = (1:size(pti,1))';
types = repmat(1, size(cam_ids,1), 1);
cams = [cam_ids, cams, types];

% object points
obj_pt_id = 1;
obj_pts = [obj_pt_id 0 0 0 2];

% image points
pti = pti + normrnd(0, 0.0001, size(pti, 1), 2);
pti_ids = (1:size(pti,1))';
cam_ids = cams(:,1);
obj_pts_ids = repmat(obj_pt_id, size(pti_ids, 1), 1);
types = repmat(1, size(pti_ids,1), 1);
pti = [pti_ids, pti, cam_ids, obj_pts_ids, types];

pte = multiray_triangulate(pti, cams, obj_pts);
norm(pte-pt)

obj_pts_res = obj_pts;
obj_pts_res(2:4) = pte;

plot_problem(1, pti, cams, obj_pts_res);
plot3(pt(:,1), pt(:,2), pt(:,3), 'g.', 'MarkerSize', 20);



