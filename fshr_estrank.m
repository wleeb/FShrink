        function khat = fshr_estrank(syf,k,gam,gam0,tol)
%
%        naive rank estimation, based on Wachter distribution threshold;
%        assumes singular values syf are sorted in decreasing order
%
        [ycut,xcut] = fshr_cuts(gam,gam0);

%%%        prin2('ycut=',ycut,1);
%%%        prin2('syf=',syf,k);
%
        khat=0;
%
        for i=1:k
%
        if (syf(i) < ycut*(1+tol))
%
        break;
    end
        khat = khat+1;
    end


        end
