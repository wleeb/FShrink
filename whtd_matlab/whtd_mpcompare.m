        function whtd_mpcompare(y,m,n,sig,ifig)
%
        s = svd(y);
        evals = s.^2;
        gam = m/n;

        whtd_mpcompare0(evals,gam,sig,ifig);

        end
%
