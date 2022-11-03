        function [fval,fgrad] = whtd_fevalvecs(tvar,what,zhat,w,z,t,m,n,k)
%
        xhat = what * tvar * zhat';
        x = w*diag(t)*z';
%
        fval = norm(xhat - x,'fro')^2 / 2;
        fgrad = what' * (xhat - x) * zhat;

        end
%
