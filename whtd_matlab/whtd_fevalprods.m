        function [fval,fgrad] = whtd_fevalprods(t_hat,t,k,...
            eout,einn,dout,dinn,cout,cinn)
%
        aa1=dout * t_hat * dinn;
        prod1 = whtd_scapro(t_hat,aa1,k*k);
%
        aa2=eout * diag(t) * einn;
        prod2 = whtd_scapro(diag(t),aa2,k*k);
%
        aa3=cout*diag(t)*cinn';
        prod3 = whtd_scapro(t_hat,aa3,k*k);
%
        fval = prod1 + prod2 - 2*prod3;
        fval = fval/2;

        fgrad = dout * t_hat * dinn - cout*diag(t)*cinn';

        end
%
