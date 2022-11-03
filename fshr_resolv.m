        function [res,der] = fshr_resolv(yeig,gam,gam0)
%
        rr = sqrt(gam+gam0-gam*gam0);
        bb = (1-gam0)*yeig - (1-gam);
        dd = bb^2 - 4*yeig*rr^2;
        res = -2 / (bb+sqrt(dd));

        bder = 1-gam0;
        dder = 2*bb*bder-4*rr^2;
        der = 2*(bder + dder/sqrt(dd)/2) / (bb+sqrt(dd))^2;


        return;

%
%        alternate formula
%

        res2 = res;
        res = ((1-gam0)*yeig + 1 - gam)^2 - 4*yeig;
        res = yeig*(1-gam0) - (1-gam) - sqrt(res);
        res = -res / 2 / (gam + gam0 - gam*gam0) / yeig;

        chk0 = res-res2;
        prin2('res=',res,1);
        prin2('res2=',res2,1);
        prin2('chk0=',chk0,1);

        prinstop;
%
        end
