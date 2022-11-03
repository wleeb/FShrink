        function [fval,fder] = mpbdry_evalf(evar,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder)
%
%        Evaluates F(lambda,e), its first partial derivative with
%        respect to e
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        fval = evar - sum(as.*awhts ./ (as*gval - rlam));

        if (ifder == 0)
%
        fder=-10000;
        return;
    end

        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

        end
%
