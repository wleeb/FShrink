        function [qval,qder] = mpbdry_evalqder(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates Q(lambda) = F(lambda,t(lambda)), where t(lambda) minimizes
%        F(lambda,e) (as a function of e); and its first derivative with
%        respect to lambda
%
        val = mpbdry_fmin(rlam,as,bs,awhts,bwhts,m,n,gam);
        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(val,rlam,as,bs,...
            awhts,bwhts,m,n,gam);

        qval = fval;
        qder = fdlam;

        end
%
