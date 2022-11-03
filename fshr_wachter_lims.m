        function [bmin,bmax] = fshr_wachter_lims(gam,gam0)
%
        rr = sqrt(gam+gam0-gam*gam0);
        bmin = (1-rr)^2 / (1-gam0)^2;
        bmax = (1+rr)^2 / (1-gam0)^2;

        end
