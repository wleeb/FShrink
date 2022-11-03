        function eleft = mpbdry_leftend(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%
%        Computes the left endpoint of the interval I_lambda; that is,
%        the root of the function G(e) - max(as) / lambda on the interval 
%        J with e > -1/(gamma * max(bs)).
%
        bmax = max(bs);
        const = min(rlam ./ as);
%
        nsteps=1000;
        vals = zeros(1,nsteps);
        es = zeros(1,nsteps);

%
%        find starting location to left of the root by bisection
%
        emin = -1 / (gam*bmax);
        einit = emin + abs(emin);

        for i=1:1000
%
        ifder=0;
        [gval,gder] = mpbdry_evalg(einit,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder);
        val = gval - const;
%
        if (val >= 0)
%
        break;
    end

        einit = (einit + emin) / 2;
    end

%
%        find root by Newton's method
%
        es(1) = einit;
%
        deps=mpbdry_machzero();
        dzero=sqrt(deps)/100;
        kbreak=0;
%
        for ijk=1:nsteps
%
        ifder=1;
        [gval,gder] = mpbdry_evalg(es(ijk),rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder);
        val = gval - const;
        der = gder;

        vals(ijk) = val;
        es(ijk+1) = es(ijk) - val / der;
%
%        check if converged
%
        if (abs(val)/const <= dzero)
%
        kbreak=kbreak+1;
    end
        if (kbreak == 2)
%
        nsteps=ijk;
        break;
    end

    end

        vals = vals(1:nsteps) / const;
        eleft = es(nsteps);


        end
%
