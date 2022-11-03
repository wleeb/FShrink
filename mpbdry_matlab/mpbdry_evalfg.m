        function [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(evar,...
            rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates F(lambda,e), its first partial derivatives with
%        respect to e; its derivative wrt lambda; and G(e) and its 
%        first derivatives
%
%
%        . . . G and its derivative (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);

%
%        derivative of F wrt lambda
%
        fdlam = -sum(as .* awhts./ (as.*gval - rlam).^2);

%
%        F and its first derivative (wrt e)
%
        fval = evar - sum(as.*awhts ./ (as*gval - rlam));
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

        end
%
