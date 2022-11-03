        function taus = fshr_taus_iid(ux,vars,m,k)
%
%        tau for each component for iid model
%
        taus = zeros(k,1);

        for j=1:k
%
        for i=1:m
%
        taus(j) = taus(j) + ux(i,j)^2 / vars(i);
    end
    end

        end
