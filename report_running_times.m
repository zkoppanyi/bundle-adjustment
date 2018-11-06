x = [25, 81, 144, 225, 256, 289, 361, 441, 529, 676]';
y = [0.748 1.093 2.7 5.76 8.114 10.501 19.257 33.462 53.69000 106.169]';

pNs = x*9; % no of unknowns


A = [x x.^2, x.^3];
coeff = A\y;
f = @(x) coeff(1).*x + coeff(2).*x.^2 + coeff(3).*x.^3;
X = linspace(min(x), max(x)+100, 100);


figure(1); clf; hold on;
plot(x, y, 'r.', 'MarkerSize', 20);
plot(X, f(X), 'r--', 'LineWidth', 2);
xlabel("No. of cameras [-]");
ylabel("Running time [s]");
set(gca, 'FontSize', 14);
grid on;

sol = fzero(@(x) f(x) - 3600/2, 1000);
fprintf("%.1f cameras take 0.5 hour\n", round(sol));

sol = fzero(@(x) f(x) - 3600, 1000);
fprintf("%.1f cameras take 1 hour\n", round(sol));

sol = fzero(@(x) f(x) - 3600*6, 1000);
fprintf("%.1f cameras take 6 hour\n", round(sol));

sol = fzero(@(x) f(x) - 3600*24, 1000);
fprintf("%.1f cameras take 24 hour\n", round(sol));