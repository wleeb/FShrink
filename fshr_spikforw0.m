        function yval = fshr_spikforw0(xval,gam,gam0)
%
        yval = (1 + xval^2) * (gam + xval^2);
        yval = yval / ((1-gam0)*xval^2 - gam0);
        yval = sqrt(yval);

        end
