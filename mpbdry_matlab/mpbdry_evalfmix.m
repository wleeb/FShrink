        function [fder2,fdmix] = mpbdry_evalfmix(evar,rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%        Evaluates the second derivative of F(lambda,e) wrt e, and its 
%        mixed partial derivative
%
%        . . . G and its first derivative (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);
%
%        second derivative of F (wrt e)
%
        term1 = sum(bs.^3 .* bwhts ./ (1+gam*evar*bs).^3);
        term2 = sum(awhts .* as.^2 ./ (as*gval - rlam).^2);
        term3 = sum(awhts .* as.^3 ./ (as*gval - rlam).^3);

        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;
%
%        mixed partial of F
%
        fdmix = 2*gder*sum(awhts .* as.^2 ./ (as*gval - rlam).^3);


        end
%
