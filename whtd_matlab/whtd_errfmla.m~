        function [err,rlam,cout,cinn] = whtd_errfmla(sx,alpha,...
            beta,gam,sig,amu,bmu)
%
%        evaluates weighted Frobenius norm error for rank 1 signal
%
        [rlam,cout,cinn] = whtd_pop2emp(sx,gam,sig);

        sout = sqrt(1 - cout^2);
        sinn = sqrt(1 - cinn^2);

        eta1 = alpha / (cout^2*alpha + sout^2*amu);
        eta2 = beta / (cinn^2*beta + sinn^2*bmu);
        eta = eta1 * eta2;

        err = sx^2*alpha*beta*(1 - cout^2*cinn^2 * eta);
        err = sqrt(err);

        end
%
