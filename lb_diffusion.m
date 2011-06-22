% LB Diffusion
% D2Q5

clear all

% constants
nw = 5;
nx = 51;
ny = 3;
niter = 1000;

tau = 1;
cs2 = 1.0/3.0;
D = cs2 * (tau - 0.5);

% lattice weights for D2Q5
w0 = 1.0/3.0;
w1 = 1.0/6.0;
w = [w0, w1, w1, w1, w1];
cx = [0, 1, 0, -1, 0];
cy = [0, 0, 1, 0, -1];

% variables
rho = zeros(1, nx, ny);

% initial conditions
fin = zeros(nw, nx, ny);
for i=1:nw
    fin(i, :, :) = w(i) * rho(1, :, :);
end

% iterate to equilibrium
for j = 1:niter
    % calculate macroscopic field
    rho = sum(fin);

    for i = 1:nw
        % collision
        feq(i, :, :) = w(i) * rho(1, :, :);
        fout(i, :, :) = fin(i, :, :) + (1.0/tau) * (feq(i, :, :) - fin(i, :, :));
        % streaming
        fin(i, :, :) = circshift(fout(i, :, :), [0, cx(i), cy(i)]);
    end

    % boundary conditions
    fin(2, 1, :) = ones(1, 1, ny) - sum(fin([1, 3, 4, 5], 1, :));
    fin(4, nx, :) = zeros(1, 1, ny) - sum(fin([1, 2, 3, 5], nx, :));
end

% Analytical comparison
delta = nx - 1;
dx = delta/(nx - 1);
x = 0:dx:delta;

rho_ana = erfc(x./(2*sqrt(D*niter)));

% plot solutions
figure(1);
hold on;
xlabel('x', 'FontSize', 16);
ylabel('rho', 'FontSize', 16);
legend('Analytical', 'LBM');
grid on;
plot(x, rho_ana, 'r o');
plot(x, rho(1,:,ny), 'b +');
