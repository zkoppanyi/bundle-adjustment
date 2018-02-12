
function show_radial_distortion_profile(prob, is_plot, tag)

    pixel_size = 3.92*10^6;
    
    if nargin < 2
        is_plot = 1;
    end
    
    if nargin < 3
        tag = '';
    end
    
    for i = 1 : size(prob.cams, 1)
        
        if size(prob.cams, 2) < 7
            fprintf('Cam #%i  f  = %.3f mm  cx = %.3f mm cy = %.3f mm\n', i, prob.cams(i, 3)*1000, prob.cams(i, 4)*1000, prob.cams(i, 5)*1000);
            continue
        end
        
        params = prob.cams(i, 6:11);
        k1 = params(1);
        k2 = params(2);
        k3 = params(3);
        p1 = params(4);
        p2 = params(5);
        
        fprintf('Cam #%i (%s) f  = %.3f mm  cx = %.3f mm cy = %.3f mm\n', i, tag, prob.cams(i, 3)*1000, prob.cams(i, 4)*1000, prob.cams(i, 5)*1000) ;
        fprintf('Cam #%i (%s) k1 = %.3e     k2 = %.3e    k3 = %.3e   p1 = %.3e  p2 = %.3e\n', i, tag, k1, k2, k3, p1, p2) ;
        
        if is_plot == 1
            r = (0 : 0.1 : 14) ; %say that's mm
            dr = k1 * r.^3 + k2 * r.^5 + k3 * r.^7;

            figure(5+i); clf; hold on;
            subplot(1,2,1); hold on;
            %plot(r, dr*1000, 'r.-');
            plot(r, dr * pixel_size / 1000, 'r.-');
            xlabel('r [mm]'); ylabel('\Deltar [\mum]');
            set(gca, 'FontSize', 12);
            grid on;
            
            ptx = -7 : 1 : 7;  % [mm]
            pty = -7 : 1 : 7;  % [mm]
            [ptx, pty] = meshgrid(ptx, pty);
            ptx = ptx(:);
            pty = pty(:);
    
            cx = 0; cy = 0;
            x_hat = (ptx - cx);
            y_hat = (pty - cy);
            r = sqrt(x_hat.^2 + y_hat.^2);
            ptix = ptx - (x_hat .* (k1 * r.^2 + k2 * r.^4 + k3 * r.^6) + p1*(r.^2 + 2*x_hat.^2) + 2*p2*x_hat.*y_hat)*10000;
            ptiy = pty - (y_hat .* (k1 * r.^2 + k2 * r.^4 + k3 * r.^6) + 2*p1*x_hat.*y_hat + p2*(r.^2 + 2*y_hat.^2))*10000;

            subplot(1,2,2); hold on;
            plot(ptx, pty, 'k.', 'MarkerSize', 7);
            plot(ptix, ptiy, 'r.', 'MarkerSize', 10);
            legend('Undistorted', 'Distorted');
            xlabel('X []'); ylabel('Y []');
            set(gca, 'FontSize', 12);
            grid on;
            axis equal;
        end

        
    end

