        function [s,der] = fshr_wachtstiel_right(x,gam,gam0)
%
        a = x*(gam+x*gam0);
        b = -1+gam+(1+gam0)*x;
%%%        c = 1;

        rad=b^2-4*a;

        s = 2 / (-b - sqrt(rad));

        da = gam + 2*gam0*x;
        db = 1+gam0;

        drad = 2*b*db - 4*da;
        bot = -b - sqrt(rad);
        dbot = -db - drad/sqrt(rad)/2;
        der = -2*dbot/bot^2;

        end
