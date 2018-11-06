function [c ceq] = const_fn(x, pts, disk_radius)
c = zeros(size(pts,1), 1);
for i = 1 : size(pts,1)
    d = (x(1) - pts(i,1)).^2 + (x(2) - pts(i,2)).^2 - x(3)^2;
    if pts(i, end) == 1
        c(i) = d;
    else
        c(i) = -d;
    end
end
ceq = [];