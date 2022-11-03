        function [cinn,cinn_der,rinn] = whtd_evalcinn(s,gam,sig)
%
%        evaluates inner cosine cinn, its derivative, and the
%        ratio between derivative and the value
%
        s0 = s/sig;
        ell0=s0^2;
        sig0=1;
        [rlam,cout,cinn] = whtd_pop2emp(s0,gam,sig0);
%
        csq = cinn^2;
        csq_der = 2*s0*(s0^4 + 2*gam*s0^2 + gam) / (s0^4 + s0^2)^2;
        cinn_der = csq_der / cinn / 2;

        cinn_der = cinn_der / sig;
%
%        ratio between derivative and value
%
        rinn = (s0^4 + 2*gam*s0^2 + gam) / (s0*(s0^2+1)*(s0^4-gam));
        rinn = rinn / sig;

        end
%
