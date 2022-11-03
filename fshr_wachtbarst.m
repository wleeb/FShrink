        function [s,sder] = fshr_wachtbarst(yeig,gam,gam0)

        [s,sder] = fshr_wachtstiel_right(yeig,gam,gam0);
        s = gam*s + gam/yeig - 1/yeig;
        sder = gam*sder - gam/yeig^2 + 1/yeig^2;


        end
