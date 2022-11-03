        function [xhat,errhat,topt] = whtd_ddexact(y,x,awhts,bwhts,m,n,k)
%
%                              description:
%
%   Evaluates optimal k-by-k diagonal matrix topt for minimizing
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
        [ux,tx,vx] = whtd_svdsmart(x,m,n,k);
        ux2 = whtd_dleft(ux,awhts,m,k);
        vx2 = whtd_dleft(vx,bwhts,n,k);
%
%        solve least squares for k-by-k topt
%
        topt = whtd_ddminvecs(uy2,vy2,tx,ux2,vx2,m,n,k,tol);
        [fval,fgrad] = whtd_ddvecs(topt,uy2,vy2,tx,ux2,vx2,m,n,k);
%
%        final matrix and error
%
        xhat = uy * diag(topt) * vy';
        errhat = sqrt(2*fval);

        end
%
