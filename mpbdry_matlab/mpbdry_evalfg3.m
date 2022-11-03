        function [fval,fder,fder2,fder3,gval,gder,gder2,gder3,fdlam,fdmix] = ...
            mpbdry_evalfg3(evar,rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates F(lambda,e), its first three partial derivatives with
%        respect to e; its derivative wrt lambda; its mixed partial; and 
%        G(e) and its first three derivatives
%
%
%        . . . G and its derivatives (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);
        gder2 = 2*gam^2*sum(bs.^3 .* bwhts./ (1 + gam*evar*bs).^3);
        gder3 = -6*gam^3*sum(bs.^4 .* bwhts./ (1 + gam*evar*bs).^4);

%
%        derivative of F wrt lambda
%
        fdlam = -sum(as .* awhts./ (as.*gval - rlam).^2);

%
%        F and its first derivative (wrt e)
%
        fval = evar - sum(as.*awhts ./ (as*gval - rlam));
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

%
%        second and third derivatives of F (wrt e)
%
        term1 = sum(bs.^3 .* bwhts ./ (1+gam*evar*bs).^3);
        dterm1 = -3*sum(bs.^4 .* bwhts ./ (1+gam*evar*bs).^4);

        term2 = sum(awhts .* as.^2 ./ (as*gval - rlam).^2);
        dterm2 = -2*gder*sum(awhts .* as.^3 ./ (as*gval - rlam).^3);


        term3 = sum(awhts .* as.^3 ./ (as*gval - rlam).^3);
        dterm3 = -3*gder*sum(awhts .* as.^4 ./ (as*gval - rlam).^4);


        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;
        fder3 = 2*gam^2*(term1*dterm2 + dterm1*term2) ...
            - 2*(2*gder*gder2*term3 + gder^2*dterm3);

%
%        mixed partial of F
%
        fdmix = 2*gder*sum(awhts .* as.^2 ./ (as*gval - rlam).^3);


        end
%
