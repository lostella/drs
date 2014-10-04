close all;
clear all;

rng(0);

% add the directory with the L1TestPack to the path
% get it from: https://www.tu-braunschweig.de/iaa/personal/lorenz/l1testpack
addpath('/media/lorenzo/Data/Datasets/l1testpack/');

m = 100; % rows of A
n = 1000; % cols of A
reg = 1;
spars = 15; % sparsity of the solution

opt.dim = [m n]; opt.sparsity = spars;
dynstr = 'dynamic2';
% run 'help construct_bpdn_instance' or 'help construct_bpdn_rhs'
% in case something doesn't work here
[A, b, xstar, tau, sigma, y] = construct_bpdn_instance('partdct', dynstr, reg, opt);
[b, y] = construct_bpdn_rhs(A, xstar, reg);
resstar = A*xstar-b;
fstar = 0.5*(resstar'*resstar)+reg*norm(xstar, 1);

stats = {};

TOL = 1e-6;
MAXIT = 1e4;

Lf = norm(A,2)^2;
mf = 0;

gammas = [0.2/Lf, (sqrt(2)-1)/Lf, 0.6/Lf, 0.8/Lf];
linest = {'-', '-.', '--', ':'};
Lh = (1+gammas*Lf)./(1-gammas*Lf).*(1-gammas*mf)./(1+gammas*mf);

%% Run DRS for different gammas

figure(1);
for i=1:length(gammas)
    optDRS.tol = TOL;
    optDRS.maxit = MAXIT;
    gam = gammas(i);
    fprintf('running DRS, gamma = %7.4e\n', gam);
    lam = 1/Lh(i);
    optDRS.fast = 0;
    t0 = tic;
    stats{end+1} = l1ls_drs(A, b, reg, gam, lam, optDRS);
    stats{end}.time = toc(t0);
    semilogy(1:stats{end}.it, (stats{end}.pobj-fstar)/(1+abs(fstar)),linest{i})
    hold on;
end

%% Run Fast DRS for different gammas

figure(2);
for i=1:length(gammas)
    optDRS.tol = TOL;
    optDRS.maxit = MAXIT;
    gam = gammas(i);
    fprintf('running Fast DRS, gamma = %7.4e\n', gam);
    lam = 1/Lh(i);
    optDRS.fast = 1;
    t0 = tic;
    stats{end+1} = l1ls_drs(A, b, reg, gam, lam, optDRS);
    stats{end}.time = toc(t0);
    semilogy(1:stats{end}.it, (stats{end}.pobj-fstar)/(1+abs(fstar)),linest{i})
    hold on;
end

