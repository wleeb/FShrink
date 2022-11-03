        function qval = mpbdry_evalq(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates Q(lambda) = F(lambda,t(lambda)), where t(lambda) minimizes
%        F(lambda,e) (as a function of e); and its first derivative with
%        respect to lambda
%
        val = mpbdry_fmin(rlam,as,bs,awhts,bwhts,m,n,gam);

        ifder=0;
        [fval,fder] = mpbdry_evalf(val,rlam,as,bs,awhts,bwhts,m,n,gam,ifder);
        qval = fval;

        end
%
