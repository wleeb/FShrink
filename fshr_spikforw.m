        function yval = fshr_spikforw(sx,gam,gam0,tau)

        xval=sqrt(tau)*sx;
%%%        xval=tau*sxf(i);
        yval = fshr_spikforw0(xval,gam,gam0);

        end
