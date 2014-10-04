close all;
clear all;

rng(0);

%% Generate problem

n    = 500; % number of variables
rc   = 1e-2; % inverse condition number of the Hessian
dens = 0.01;

% Random QP generator taken from AN OPTIMAL ALGORITHM FOR CONSTRAINED DIFFERENTIABLE CONVEX OPTIMIZATION
% GONZAGA, ELIZABETH W. KARAS, AND DIANE R. ROSSETTO
Q  = sprandsym(n,dens,rc,1);% Hessian
xmin = -1e2*ones(n,1);% lower bound
xmax = 1e1*ones(n,1);% upper bound
z = 1e2*randn(n,1);
xstar = min(max(z,xmin),xmax);
q2 = zeros(n,1);
q2(xstar==xmin) = rand(sum(xstar==xmin),1);
q2(xstar==xmax) = -rand(sum(xstar==xmax),1);
q2 = 1e1*q2;
q = q2-Q*xstar;
fstar = 1/2*xstar'*Q*xstar+q'*xstar;

stats = {};

TOL = 1e-6;
MAXIT = 1e4;

eigQ = eig(Q);
Lf = max(eigQ);
mf = min(eigQ);

gammas = [0.2/Lf, (sqrt(2)-1)/Lf, 0.6/Lf, 0.8/Lf];
linest = {'-', '-.', '--', ':'};
Lh = (1+gammas*Lf)./(1-gammas*Lf).*(1-gammas*mf)./(1+gammas*mf);

%% QUADPROG
%optquadprog.Algorithm = 'interior-point-convex';
%optquadprog.Display = 'none';
%t0 = tic;
%[x,f,exitflag,output] = quadprog(Q, q, [], [], [], [], xmin, xmax, [], optquadprog);
%stats{end+1}.time = toc(t0);
%stats{end}.pobj = f;
%stats{end}.x = x;

%% Run DRS
figure(1);
for i=1:length(gammas)
    optDRS.tolr = TOL;
    optDRS.maxit = MAXIT;
    optDRS.fast = 0;
    gam = gammas(i);
    lam = 1/Lh(i);
    fprintf('Running DRS, gamma = %7.4e\n', gam);
    t0 = tic;
    stats{end+1} = boxqp_drs(Q, q, xmin, xmax, gam, lam, optDRS);
    stats{end}.time = toc(t0);
    semilogy(1:stats{end}.it, (stats{end}.pobj-fstar)/(1+abs(fstar)),linest{i})
    hold on;
end

%% Run DRS
figure(2);
for i=1:length(gammas)
    optDRS.tolr = TOL;
    optDRS.maxit = MAXIT;
    optDRS.fast = 1;
    gam = gammas(i);
    lam = 1/Lh(i);
    fprintf('Running Fast DRS, gamma = %7.4e\n', gam);
    t0 = tic;
    stats{end+1} = boxqp_drs(Q, q, xmin, xmax, gam, lam, optDRS);
    stats{end}.time = toc(t0);
    semilogy(1:stats{end}.it, (stats{end}.pobj-fstar)/(1+abs(fstar)),linest{i})
    hold on;
end
