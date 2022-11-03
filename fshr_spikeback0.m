        function xval = fshr_spikeback0(yval,gam,gam0)
%
        bb = 1 + gam - yval^2*(1-gam0);
        cc = gam + yval^2*gam0;

        xval = -bb + sqrt(bb^2 - 4*cc);
        xval = xval / 2;
        xval = sqrt(xval);


        return;


        [val,der] = fshr_evalpsi(yval^2,gam,gam0);
        xval2 = 1/sqrt(val);

        chk0 = abs(xval2 - xval);
        prin2('chk0=',chk0,1);

        prinstop;

        end
