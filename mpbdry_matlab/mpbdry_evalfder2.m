        function [fder,fder2] = mpbdry_evalfder2(evar,rlam,as,bs,awhts,...
            bwhts,m,n,gam,ifder2)
%
%        Evaluates the first and second partial derivatives of F(lambda,e) with
%        respect to e; second derivative is returned only if ifder2 is not 0
%
%
%        . . . G and its derivative (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);

%
%        first derivative of F (wrt e)
%
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

        if (ifder2 == 0)
%
        fder2=-1000;
        return
    end
%
%        second derivative of F (wrt e)
%
        term1 = sum(bs.^3 .* bwhts ./ (1+gam*evar*bs).^3);
        term2 = sum(awhts .* as.^2 ./ (as*gval - rlam).^2);
        term3 = sum(awhts .* as.^3 ./ (as*gval - rlam).^3);

        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;

        end
%
