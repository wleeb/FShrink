        function xcut = fshr_xcut_tau(gam,gam0,tau)

        [ycut,xcut] = fshr_cuts(gam,gam0);
        xcut = xcut / sqrt(tau);

        end
