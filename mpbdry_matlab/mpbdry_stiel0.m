        function [s,sder,evar,eder,gval,gder,gderl] = mpbdry_stiel0(rlam,...
            as,bs,awhts,bwhts,m,n,gam)

        x0=0;
        [xs,vals,nsteps] = mpbdry_rootf(x0,rlam,as,bs,awhts,bwhts,m,n,gam);
        evar = xs(nsteps);

        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(evar,rlam,as,bs,awhts,...
            bwhts,m,n,gam);

        eder = -fdlam / fder;
        gderl = gder * eder;

%
%        evaluate the stieltjes transform and derivative
%
        s = sum(1 ./ (gval*as-rlam).*awhts);
        sder = -sum(awhts.*(gderl*as - 1) ./ (gval*as - rlam).^2);


        return

%
%        evaluate them the slow way
%
        s = 0;
        for i=1:m
%
        s = s + 1 / (as(i)*gval - rlam) * awhts(i);
    end

        sder = 0;
        for i=1:m
%
        sder = sder - (as(i)*gderl - 1) / (as(i)*gval - rlam)^2 * awhts(i);
    end

        end
%
