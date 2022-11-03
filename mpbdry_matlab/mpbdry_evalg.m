        function [gval,gder] = mpbdry_evalg(evar,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder)
%
%        Evaluates G(e) and its first derivative
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));

        if(ifder==0)
%
        gder=-10000;
        return;
    end

        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);

        end
%
