        function [fval,fgrad] = whtd_ddprods(tvec,tx,dout,dinn,cout,cinn,...
            eout,einn,m,n,k)
%
        fgrad = zeros(k,1);

        prin2('tvec=',tvec,k);
        prin2('tx=',tx,k);

        aa1=dout * diag(tvec) * dinn;
        prod1 = whtd_scapro(diag(tvec),aa1,k*k);
%
        aa2=eout * diag(tx) * einn;
        prod2 = whtd_scapro(diag(tx),aa2,k*k);
%
        aa3=cout*diag(tx)*cinn';
        prod3 = whtd_scapro(diag(tvec),aa3,k*k);
%
        fval = prod1 + prod2 - 2*prod3;
        fval = fval/2;


        aa = dout.*dinn;
        bb = cout.*cinn;
        fgrad = aa*tvec - bb*tx;

        end
%
