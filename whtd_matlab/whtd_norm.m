        function wnorm = whtd_norm(x,awhts,bwhts,m,n)
%
        x2 = whtd_dmult(x,awhts,bwhts,m,n);
        wnorm = norm(x2,'fro');

        end
%
