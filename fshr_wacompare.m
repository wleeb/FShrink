        function fshr_wacompare(yf,gam,gam0,ifig)
%
%        plots histogram of squared singular values of yf 
%        and overlays plot of continuous Wachter density for white noise 
%        with variance 1, with aspect ratios gam and gam0.
%
        s_all = svd(yf);
        evals = s_all.^2;
        fshr_wacompare0(evals,gam,gam0,ifig);

        end
