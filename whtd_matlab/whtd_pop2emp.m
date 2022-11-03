        function [sy,cout,cinn] = whtd_pop2emp(tx,gam,sig)
%
        tx=tx/sig;
%
        tmin = sqrt(sqrt(gam));
        if (tx <= tmin);
%
        sy=sig*(1 + sqrt(gam));
        cout=0;
        cinn=0;
        return;
    end

        sy = (tx^2 + 1)*(1 + gam/tx^2);
        sy = sig*sqrt(sy);

        cout = (tx^4 - gam) / (tx^4 + tx^2*gam);
        cout = sqrt(cout);

        cinn = (tx^4 - gam) / (tx^4 + tx^2);
        cinn = sqrt(cinn);

        end
%
