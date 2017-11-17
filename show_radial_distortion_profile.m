
function show_radial_distortion_profile(prob, is_plot)

    if nargin == 1
        is_plot = 1;
    end
    
    for i = 1 : size(prob.cams, 1)
        
        params = prob.cams(i, 6:11);
        dc = params(1);
        k1 = params(2);
        k2 = params(3);
        k3 = params(4);
        p1 = params(5);
        p2 = params(6);

        if is_plot == 1,
            r = 0 : 0.1 : 14;
            dr = k1 * r.^3 + k2 * r.^5 + k3 * r.^7;

            figure(5+i); clf; hold on;
            plot(r, dr*1000, 'r.-');
            xlabel('r [mm]'); ylabel('\Deltar [\mum]');
            set(gca, 'FontSize', 12);
            grid on;
        end

        fprintf("Cam #%i  %.3f %.3e %.3e %.3e %.3e %.3e\n", i,  dc, k1, k2, k3, p1, p2) ;
    end

