        function f = fshr_evalwacht(x,gam,gam0)
%
        [bmin,bmax] = fshr_wachter_lims(gam,gam0);

        dd = (1-gam0) * sqrt((bmax - x) * (x - bmin));
        dn = 2*pi*x*(gam + x*gam0);
        f = dd/dn;

        end
