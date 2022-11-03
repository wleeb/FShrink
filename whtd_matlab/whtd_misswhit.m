        function [yw,awhts,bwhts] = whtd_misswhit(yholes,qr,qc,m,n)
%
        awhts = sqrt(1./qr);
        bwhts = sqrt(1./qc);
%
%        equalize the noise variances
%
        yw = whtd_dmult(yholes,awhts,bwhts,m,n);

        end
%
