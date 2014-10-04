function result = boxqp_drs(Q, q, lx, ux, gam, lam, opt)

n = size(Q,1);

if nargin < 7 || ~isfield(opt, 'maxit')
    maxit = 1e4;
else
    maxit = opt.maxit;
end
if nargin < 7 || ~isfield(opt, 'tolr')
    tolr = 1e-4;
else
    tolr = opt.tolr;
end
if nargin < 7 || ~isfield(opt, 'ptarget')
    checkpobj = 0;
else
    checkpobj = 1;
    ptarget = opt.ptarget;
end
if nargin < 7 || ~isfield(opt, 'fast')
    fast = 0;
else
    fast = opt.fast;
end
if nargin < 7 || ~isfield(opt, 'x0')
    x = (lx+ux)/2;
else
    x = opt.x0;
end

result.ts = zeros(1,maxit);
result.pobj = zeros(1,maxit);
result.pres = zeros(1,maxit);

t0 = tic;

L = chol(eye(n)+gam*Q,'lower');
Lt = L';

tau = 1;
x1 = x;
w = x;

for it=1:maxit
    % y-step
    y = Lt\(L\(w-gam*q));
    % z-step
    z = max(lx,min(ux,2*y-w));
    % x-step
    x = w+lam*(z-y);
    
    pobj = 0.5*(z'*(Q*z)) + q'*z;
    
    % record iterates
    result.ts(1,it) = toc(t0);
    result.pobj(1,it) = pobj;
    result.pres(1,it) = norm(z-y);
    
    % stopping criterion
    if checkpobj == 0 && result.pres(1,it) <= tolr
        break;
    elseif checkpobj == 1 && pobj <= ptarget
        break;
    end
    
    % extrapolation
    if fast
        beta = it/(it+3);
    else
        beta = 0;
    end
    w = x+beta*(x-x1);
    x1 = x;
end
result.it = it;
result.ts = result.ts(1,1:it);
result.pobj = result.pobj(1,1:it);
result.pres = result.pres(1,1:it);
result.x = x;
result.y = y;
result.z = z;
