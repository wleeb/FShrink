        function [val,der] = fshr_evalpsi(yeig,gam,gam0)
%
        [s,sder] = fshr_wachtbarst(yeig,gam,gam0);
        [r,rder] = fshr_resolv(yeig,gam,gam0);
%
        val= s*r*yeig;
        der = s*r + yeig*sder*r + yeig*s*rder;


        return

%
%        alternative formula for psi
%
        val2=val;

        bb = (1-gam0)*yeig - 1 - gam;
        ff = (1-gam0)*yeig + 1 - gam;
        dd = ff^2  - 4*yeig;
        aa = 2*(gam0*yeig + gam);
%
        val = (bb - sqrt(dd)) / aa;

        chk0 = val2 - val;
        prin2('chk0=',chk0,1);

        prinstop;

        end
