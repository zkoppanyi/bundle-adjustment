function show_sparisty_pattern(J, fig_num)

    sparI = J;
    sparI(sparI < 1e-4) = 0;
    sparI(sparI == 0) = 0;
    sparI(sparI ~= 0) = 255;

    figure(fig_num); clf; hold on;
    image(flip(sparI));
    grid on;
    axis equal;