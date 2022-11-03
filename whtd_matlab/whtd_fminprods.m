        function tmin = whtd_fminprods(t,k,dout,dinn,cout,cinn,tol)
%
        [dout_inv,khat] = whtd_pseudoinv(dout,k,k,tol);
        [dinn_inv,khat] = whtd_pseudoinv(dinn,k,k,tol);
        tmin = dout_inv*cout*diag(t)*cinn'*dinn_inv';

        end
%
