        function [cout,cout_der,rout] = whtd_evalcout(s,gam,sig)
%
%        evaluates outer cosine cout, its derivative, and the
%        ratio between derivative and the value
%
        s0 = s/sig;
        ell0=s0^2;
        sig0=1;
        [rlam,cout,cinn] = whtd_pop2emp(s0,gam,sig0);
%
        csq = cout^2;
        csq_der = 2*gam*s0*(s0^4 + 2*s0^2 + gam) / (s0^4 + gam*s0^2)^2;
        cout_der = csq_der / cout / 2;
        cout_der = cout_der / sig;
%
%        ratio between derivative and value
%
        rout = gam*(s0^4 + 2*s0^2 + gam) / (s0 * (s0^2+gam) * (s0^4-gam));
        rout = rout / sig;

        end
%
