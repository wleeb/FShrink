        function uptwo = fshr_uptwo_fmla(z,gam,gam0)
%
        [res,dres] = fshr_resolv(z,gam,gam0);
        [sbar,sbder] = fshr_wachtbarst(z,gam,gam0);

        [s,sder] = fshr_mpstiel_left(-sbar,gam0);

        ff = 1 + gam*res + gam*z*dres;
        ff = ff / z^2;

        gg = 1 - 2*sbar*s + sbar^2*sder;

        uptwo = ff*gg;




%%%        prinstop;

        end
