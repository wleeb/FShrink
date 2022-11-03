        function tau = fshr_esttau_iid(ep0,m,n0)
%
        gam0 = m/n0;

        sep0 = svd(ep0);
        tau = sum(1./sep0.^2) / m  * (1-gam0);

%%%        tau = 1 / tau / (1-gam0);
%%%        prin2('tau, estimated=',tau,1);


        end
