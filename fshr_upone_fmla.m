        function [upone,dupone] = fshr_upone_fmla(z,gam,gam0)
%

        [sbar,sbder] = fshr_wachtbarst(z,gam,gam0);
        [ss,ssder] = fshr_mpstiel_left(-sbar,gam0);


        upone = 1/z - sbar/z + sbar^2*ss/z;

        gg = sbar^2*ss;
        ff = 1 - sbar + gg;
%%%        upone = ff/z;

        dss=-ssder*sbder;
        dgg=2*sbar*sbder*ss + sbar^2*dss;
        dff = -sbder + dgg;
        dupone = dff/z - ff/z^2;

        end
