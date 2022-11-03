        function [ycut,xcut] = fshr_cuts(gam,gam0)
%
%        singular value thresholds
%
        ycut = 1 + sqrt(gam0 + gam - gam*gam0);
        ycut = ycut / (1-gam0);

        xcut = gam0 + sqrt(gam0 + gam - gam*gam0);
        xcut = xcut / (1-gam0);
        xcut = sqrt(xcut);


        end
