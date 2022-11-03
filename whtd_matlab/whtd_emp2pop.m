        function [tx,cout,cinn] = whtd_emp2pop(rlam,gam,sig)
%
        rlam=rlam/sig;

        btop = 1+sqrt(gam);
        if (rlam <= btop);
%
        tx=0;
        cout=0;
        cinn=0;
        return;
    end

        yy = rlam^2 - 1 - gam;
        tx = (yy + sqrt(yy^2 - 4*gam))/2;
        tx = sqrt(tx);
%
        cout = (tx^4 - gam) / (tx^4 + tx^2*gam);
        cout = sqrt(cout);
%
        cinn = (tx^4 - gam) / (tx^4 + tx^2);
        cinn = sqrt(cinn);

        tx=tx*sig;

        end
%
