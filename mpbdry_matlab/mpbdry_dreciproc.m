        function [val,vder] = mpbdry_dreciproc(rlam,as,bs,awhts,bwhts,m,n,gam)
%
        [s,sder] = mpbdry_stiel(rlam,as,bs,awhts,bwhts,m,n,gam);
        [sbar,sbar_der] = mpbdry_stra2sbar(s,sder,rlam,gam);
        d = s * sbar * rlam;
%
%        . . . derivative of D wrt rlam, NOT sqrt(rlam)
%
        dder = sder*sbar*rlam + s*sbar_der*rlam + s*sbar;
%
%        1/D, and its derivative
%
        val = 1/d;
        vder = -dder/d^2;

        end
%
