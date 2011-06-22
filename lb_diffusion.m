% LB Poiseuille Flow
% D2Q9

clear all

% constants
global nw;
nw = 9;
nx = 31;
ny = 3;
niter = 4000;

global tau;
global cs2;
global F;
tau = 1;
cs2 = 1.0/3.0;
Rho = 1;
F = [0, 0.0002];
nu = cs2 * (tau - 0.5);
mu = nu * Rho;

% lattice weights
w0 = 4.0/9.0;
w1 = 1.0/9.0;
w2 = 1.0/36.0;
global w;
global cx;
global cy;
w = [w0, w1, w1, w1, w1, w2, w2, w2, w2];
cx = [0, 1, 0, -1, 0, 1, -1, -1, 1];
cy = [0, 0, 1, 0, -1, 1, 1, -1, -1];
oppositeOf = [1, 4, 5, 2, 3, 8, 9, 6, 7];

% variables
global rho;
global ux;
global uy;
global fin;

% initial conditions
fin = zeros(nw, nx, ny);
for i=1:nw
    fin(i,:,:) = pois_feq(i, ones(Rho, nx, ny), 0, 0);
end

% iterate to equilibrium
for i=1:niter
    % calculate macroscopic field
    [rho, ux, uy] = my_rho(rho, ux, uy);

    % collision
    for i=1:nw
        feq(i, :, :) = pois_feq(i, rho, ux, uy);
        fout(i,:,:) = fin(i,:,:) + (1/tau) * (feq(i, :, :) - fin(i,:,:)) + (w(i)/cs2)*(cx(i)*F(1) + cy(i)*F(2));

        % bounce back
        fout(i,[1 nx],:) = fin(oppositeOf(i),[1 nx],:);

        % streaming
        fin(i,:,:) = circshift(fout(i,:,:), [0 cx(i) cy(i)]);
    end

end

% Analytical comparison
figure(1);
delta = nx - 1;
dx = delta/(nx - 1);
x = 0:dx:delta;
uy_ana = -(0.5/mu)*F(2)*(x.*x - delta.*x + (delta-dx/2.0)*dx/2);
hold on;
plot(x, uy_ana, 'r o');
plot(x, uy(1,:,ny), 'b +');
xlabel('x', 'FontSize', 16);
ylabel('uy', 'FontSize', 16);
legend('Analytical', 'LBM');
grid on;
