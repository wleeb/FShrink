        function [qval,qder,qder2] = mpbdry_evalqder2(rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%        Evaluates Q(lambda) = F(lambda,t(lambda)), where t(lambda) minimizes
%        F(lambda,e) (as a function of e); and its first two derivatives with
%        respect to lambda
%
        [val,der] = mpbdry_fminder(rlam,as,bs,awhts,bwhts,m,n,gam);
        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(val,rlam,as,bs,...
            awhts,bwhts,m,n,gam);

        qval = fval;
        qder = fdlam;
%
        qder2 = 2*(as.*awhts.*(as*gder*der-1) ./ (gval*as - rlam).^3);
        qder2=sum(qder2);

        return
%
%        compute qder2 the slow way:
%
        qder2=0;
        for i=1:m
%
        qder2 = qder2 + 2*as(i)*(as(i)*gder*der - 1) / ...
            (as(i)*gval - rlam)^3 * awhts(i);
    end

        end
%
