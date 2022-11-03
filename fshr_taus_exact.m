        function taus = fshr_taus_exact(dd,vars,m,k)
%
        taus = zeros(k,1);

        for j=1:k
%
        for i=1:m
%
        taus(j) = taus(j) + dd(i,j)^2 / vars(i) / m;
    end
    end

%%%        prin2('taus=',taus,k);

        end
