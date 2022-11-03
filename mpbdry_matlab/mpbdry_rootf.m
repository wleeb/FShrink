        function [xs,fvals,nsteps] = mpbdry_rootf(x0,rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%        Finds the root of F_lambda(e), as a function of e, by Newton
%
        nsteps = 1000;
        xs = zeros(1,nsteps);
        fvals = zeros(1,nsteps);
        fders = zeros(1,nsteps);
%
%        initialize at x0
%
        xs(1)=x0;

        deps = mpbdry_machzero();
%%%        tol=100*deps;
        tol=sqrt(deps)/10;

        kstop=0;
        for i=1:nsteps
%
        ifder=1;
        [fvals(i),fders(i)] = mpbdry_evalf(xs(i),rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder);
        xs(i+1) = xs(i) - fvals(i) / fders(i);

        if (abs(fvals(i)) < tol)
%
        kstop = kstop+1;
    end
        if (kstop == 2)
%
        nsteps = i;
        break;
    end
    end
%
        fvals = fvals(1:nsteps);
        xs = xs(1:nsteps);

        end
%
