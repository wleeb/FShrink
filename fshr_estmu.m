        function dmu = fshr_estmu(cov0,m)

        dmu = 0;
        for i=1:m
%
        dmu = dmu + cov0(i,i) / m;
    end

        end
