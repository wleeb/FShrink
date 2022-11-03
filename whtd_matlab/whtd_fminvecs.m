        function tmin = whtd_fminvecs(what,zhat,w,z,t,m,n,k,tol)
%
        [winv,khat] = whtd_pseudoinv(what,m,k,tol);
        [zinv,khat] = whtd_pseudoinv(zhat,n,k,tol);
        tmin = winv*w*diag(t)*z'*zinv';

        end
%
