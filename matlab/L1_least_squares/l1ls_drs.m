function result = l1ls_drs(A, b, reg, gam, lam, opt)

    m = size(A,1);
    n = size(A,2);

    if nargin < 6 || ~isfield(opt, 'maxit')
        maxit = 1e4;
    else
        maxit = opt.maxit;
    end
    if nargin < 6 || ~isfield(opt, 'tolr')
        tolr = 1e-4;
    else
        tolr = opt.tol;
    end
    if nargin < 6 || ~isfield(opt, 'ptarget')
        checkpobj = 0;
    else
        checkpobj = 1;
        ptarget = opt.ptarget;
    end
    if nargin < 6 || ~isfield(opt, 'fast')
        fast = 0;
    else
        fast = opt.fast;
    end
    if nargin < 6 || ~isfield(opt, 'x0')
        x = zeros(n, 1);
    else
        x = opt.x0;
    end

    In = eye(n);
    result.ts = zeros(1,maxit);
    result.pobj = zeros(1,maxit);
    result.pres = zeros(1,maxit);

    t0 = tic;

    L = chol(A'*A+1/gam*In,'lower');
    Lt = L';
    Atb = A'*b;

    tau = 1;
    x1 = x;
    w = x;

    for it=1:maxit
        % y-step
        y = Lt\(L\(Atb+w/gam));
        % z-step
        z1 = 2*y-w;
        z = sign(z1).*max(0,abs(z1)-gam*reg);
        % x-step
        x = w+lam*(z-y);
        
        res = A*z-b;
        pobj = 0.5*(res'*res) + reg*norm(z,1);
        
        % record iterates
        result.ts(1,it) = toc(t0);
        result.pobj(1,it) = pobj;
        result.pres(1,it) = norm(z-y);
        
        % stopping criterion
        if checkpobj == 0 && result.pres(1,it) < tolr
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
    result.z = z;
end
