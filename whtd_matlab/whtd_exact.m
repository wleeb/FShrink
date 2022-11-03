        function [xhat,errhat,topt] = whtd_exact(y,x,awhts,bwhts,m,n,k)
%
%                              description:
%
%   Evaluates optimal k-by-k matrix topt for minimizing
%
%                || A (X - \hat U topt \hat V^T) B^T ||_{Fr}^2      (1)
%
%   where A and B are the diagonal weight matrices, and \hat U and \hat V
%   are the empirical singular vectors of Y.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   x - the m-by-n matrix of data vectors
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - the Frobenius loss between X and \hat X
%   topt - optimal k-by-k matrix used in spectral denoiser
%   
%
%        . . . apply weight matrices to empirical vectors
%

        tol = 100*eps;
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);

%
%        apply weight matrices to population vectors
%
        [ux,sx,vx] = whtd_svdsmart(x,m,n,k);
        ux2 = whtd_dleft(ux,awhts,m,k);
        vx2 = whtd_dleft(vx,bwhts,n,k);

%
%        solve least squares for k-by-k topt
%
        topt = whtd_fminvecs(uy2,vy2,ux2,vx2,sx,m,n,k,tol);
        [fval,fgrad] = whtd_fevalvecs(topt,uy2,vy2,ux2,vx2,sx,m,n,k);
%%%        chk0 = norm(uy2 * topt * vy2' - ux2 * diag(sx) * vx2','fro') ...
%%%            - sqrt(2*fval)
%
%        final matrix and error
%
        xhat = uy * topt * vy';
        errhat = sqrt(2*fval);
%%%        err = norm(diag(awhts)*(xhat - x)*diag(bwhts),'fro');
%%%        chk0 = errhat - err

        end
%
