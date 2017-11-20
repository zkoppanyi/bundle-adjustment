clear variables; clc; close all;

ptx = -7360/2 : 300 : 7360/2;
pty = -4912/2: 300 : 4912/2;
[ptx, pty] = meshgrid(ptx, pty)
ptx = ptx(:) * 4.87e-3;
pty = pty(:) * 4.87e-3;

cx = 0;
cy = 0;

cam.k1 = 5e-6;
cam.k2 = 5e-8;
cam.k3 = 5e-10;
cam.p1 = 5e-6;
cam.p2 = 5e-6;

x_hat = ptx - cx;
y_hat = pty - cy;
r = sqrt(x_hat.^2 + y_hat.^2);
ptix = ptx - x_hat .* (cam.k1 * r.^2 - cam.k2 * r.^4 - cam.k3 * r.^6) - cam.p1*(r.^2 - 2*x_hat.^2) - 2*cam.p2*x_hat.*y_hat;
ptiy = pty - y_hat .* (cam.k1 * r.^2 - cam.k2 * r.^4 - cam.k3 * r.^6) - 2*cam.p1*x_hat.*y_hat - cam.p2*(r.^2 - 2*y_hat.^2);
        
figure(1); hold on;
plot(ptx, pty, 'r.', 'MarkerSize', 10);
plot(ptix, ptiy, 'g.', 'MarkerSize', 10);

img_pts = [ptix, ptiy];
img_pts_ud = [ptx, pty];
dr = img_pts(:,1:2) - img_pts_ud(:,1:2);
dr = sqrt(dr(:,1).^2 + dr(:,2).^2);
r = sqrt(img_pts(:,1).^2 + img_pts(:,2).^2);
figure(5); clf; hold on;
plot(r, dr, 'r*')

    
