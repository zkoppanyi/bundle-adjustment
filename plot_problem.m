
function plot_problem(plot_k, img_pts, imgs, obj_pts, cams, color_in)

    %figure('WindowScrollWheelFcn',@figScroll);
    %clf; 
    hold on;

    color = 'r';
    if nargin > 5
        color = color_in;
    end
    
    %cam_scale = 1;
    cam_scale = 20;
    
    % plot cameras
    Rcams = {};
    for i = 1 : size(imgs, 1)
        
        cam_id = imgs(i, 8);
        f  = cams(cam_id,3);
        cx = cams(cam_id,4);
        cy = cams(cam_id,5);

        cam_x  = imgs(i, 2);
        cam_y  = imgs(i, 3);
        cam_z  = imgs(i, 4);

        omega = imgs(i, 5);
        phi   = imgs(i, 6);
        kappa = imgs(i, 7);

        plot3(cam_x, cam_y, cam_z, 'r*');
        R = get_rotation_matrix(omega, phi, kappa);
        cam_dir = [0 0 f] * R * cam_scale;
        plot3([cam_x cam_x+cam_dir(1)], [cam_y cam_y+cam_dir(2)], [cam_z cam_z+cam_dir(3)], [color '-']);

        % camera plot
        pl = [-1 -1 0; -1 1 0; 1 1 0; 1 -1 0; -1 -1 0];
        cam_dir_f = [0 0 f] * R * cam_scale;
        plt = pl*R * cam_scale + repmat([cam_x - cam_dir_f(1), cam_y - cam_dir_f(2), cam_z - cam_dir_f(3)], size(pl,1), 1);
        plot3(plt(:,1), plt(:,2), plt(:,3), [color '.-']);
        for k = 1 : size(plt,1)
            plot3([cam_x plt(k,1)], [cam_y plt(k,2)], [cam_z plt(k,3)], [color '.-']);
        end
        
        Rcams{i} = R;
    end
    
    for i = 1 : size(img_pts, 1)

        img_id = img_pts(i, 4);
        cam_id = imgs(img_id, 8);
        obj_pt_id = img_pts(i, 5);     
        
        f  = cams(cam_id,3);
        cx = cams(cam_id,4);
        cy = cams(cam_id,5);
        cam_x  = imgs(img_id, 2);
        cam_y  = imgs(img_id, 3);
        cam_z  = imgs(img_id, 4);
        R = Rcams{img_id};
 
        plot3(obj_pts(obj_pt_id,2), obj_pts(obj_pt_id,3), obj_pts(obj_pt_id,4), 'b.', 'MarkerSize', 20);
        plot3([cam_x obj_pts(obj_pt_id,2)], [cam_y obj_pts(obj_pt_id,3)], [cam_z obj_pts(obj_pt_id,4)], 'b-');

        pti2 = [img_pts(i,2)-cx, img_pts(i,3)-cy, -f] * R;
        plot3(cam_x + pti2(1) * cam_scale, cam_y + pti2(2) * cam_scale, cam_z + pti2(3) * cam_scale, 'r.','MarkerSize', 12);
    end
    
    xlabel("X"); ylabel("Y"); zlabel("Z");
    grid on;
    axis equal;
    
    function figScroll(src, callbackdata)
        disp("scroll")
%       if callbackdata.VerticalScrollCount > 0 
%          xd = h.XData;
%          % This code uses dot notation to set properties
%          % Dot notation runs in R2014b and later.
%          % For R2014a and earlier: xd = get(h,'XData');
%          inc = xd(end)/20;
%          x = [0:.1:xd(end)+inc];
%          re_eval(x)
%       elseif callbackdata.VerticalScrollCount < 0 
%          xd = h.XData;
%          % For R2014a and earlier: xd = get(h,'XData');
%          inc = xd(end)/20;
%          x = [0:.1:xd(end)-inc+.1]; % Don't let xd = 0;
%          re_eval(x)
%       end
    end

end
