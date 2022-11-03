        function [s,der] = fshr_mpstiel_left(z,gam)
%
%        MP stieltjes transform, valid for z on left of support
%
        a = gam*z;
        b = z + gam - 1;
        s = 2 / (-b + sqrt(b^2 - 4*a));

        rad = sqrt((z - 1 -gam)^2 -4*gam);
        drad = 2*(z - 1 -gam)/rad/2;

        bot = 1-gam-z + rad;
        dbot = -1 + drad;
        der = -2*dbot/bot^2;


        end
