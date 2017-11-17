

for i = 1 : size(prob.cams, 1)
    params = sol.cams(i, 6:10);
    k1 = params(1);
    k2 = params(2);
    k3 = params(3);
    p1 = params(4);
    p2 = params(5);
    
    r = 0 : 0.1 : 14;
    dr = k1 * r.^3 + k2 * r.^5 + k3 * r.^7;
    
    figure(5+i); clf; hold on;
    plot(r, dr*1000, 'r.-');
    xlabel('r [mm]'); ylabel('\Deltar [\mum]');
    set(gca, 'FontSize', 12);
    grid on;
    
    fprintf("Cam #%i  %.3e %.3e %.3e %.3e %.3e\n", i,  params(1), params(2), params(3), params(4), params(5)) ;
end

