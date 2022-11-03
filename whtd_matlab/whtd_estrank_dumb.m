        function kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin)
%
        rlams = sy.^2;
%
        bedge = sig^2*(1 + sqrt(gam))^2;
        bedge = bedge + sig^2*bmargin;
        kest=0;
        for i=1:k
%
        if (rlams(i) <= bedge)
            break;
        end
        kest = kest+1;
    end

        end
%
