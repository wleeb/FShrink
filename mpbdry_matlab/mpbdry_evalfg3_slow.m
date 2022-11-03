        function [fval,fder,fder2,fder3,gval,gder,gder2,gder3,fdlam,fdmix] = ...
            mpbdry_evalfg3_slow(evar,rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates F(lambda,e), its first three partial derivatives with
%        respect to e; its derivative wrt lambda; its mixed partial; and 
%        G(e) and its first three derivatives
%
%
%        . . . G and its derivatives (wrt e)
%
        gval = 0;
        gder=0;
        gder2 =0;
        gder3=0;
        for j=1:n
%
        gval = gval + bs(j) / (1 + gam*bs(j)*evar) * bwhts(j);
        gder = gder - gam*bs(j)^2 / (1 + gam*bs(j)*evar)^2  * bwhts(j);
        gder2 = gder2 + 2*gam^2*bs(j)^3 / (1 + gam*bs(j)*evar)^3  * bwhts(j);
        gder3 = gder3 - 6*gam^3*bs(j)^4 / (1 + gam*evar*bs(j))^4 * bwhts(j);
    end

%
%        derivative of F wrt lambda
%
        fdlam = 0;
        for j=1:m
%
        fdlam = fdlam - as(j) / (as(j)*gval - rlam)^2 * awhts(j);
    end

%
%        F and its first derivative (wrt e)
%
        fval = 0;
        fder = 0;
        for i=1:m
%
        fval = fval + as(i) / (as(i)*gval - rlam) * awhts(i);
        fder = fder - gder*as(i)^2 / (as(i)*gval - rlam)^2  * awhts(i);
    end
        fval = evar - fval;
        fder = 1-fder;

%
%        second and third derivatives of F (wrt e)
%
        term1 = 0;
        term2 = 0;
        term3 = 0;

        dterm1=0;
        dterm2=0;
        dterm3=0;

        for i=1:n
%
        term1 = term1+bs(i)^3/(1 + gam*bs(i)*evar)^3 * bwhts(i);
        dterm1 = dterm1-3*gam*bs(i)^4/(1 + gam*bs(i)*evar)^4 * bwhts(i);
    end
%
        for i=1:m
%
        term2 = term2 + as(i)^2 / (as(i) * gval - rlam)^2 * awhts(i);
        dterm2 = dterm2 -2*gder* as(i)^3 / (as(i) * gval - rlam)^3 * awhts(i);

        term3 = term3 + as(i)^3 / (as(i)*gval - rlam)^3 * awhts(i);
        dterm3 = dterm3 - 3*gder*as(i)^4 / (as(i)*gval - rlam)^4 * awhts(i);
    end
%
        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;
        fder3 = 2*gam^2*(term1*dterm2 + dterm1*term2) ...
            - 2*(2*gder*gder2*term3 + gder^2*dterm3);

%
%        mixed partial of F
%
        fdmix=0;
        for i=1:m
%
        fdmix = fdmix + 2*gder*as(i)^2 / (as(i)*gval - rlam)^3 * awhts(i);
    end


        end
%
