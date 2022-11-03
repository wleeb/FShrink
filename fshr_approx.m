        function [xhat,xfhat,errhat,khat,uhat,vhat,tvals,etas,taus,sx,couts,...
            cinns,dmu] = fshr_approx(y,cov0,gam0,m,n,k,iorth)
%
%                                description:
%
%   Performs asymptotically optimal singular value shrinkage with pseudo-whitening,
%   for input of the form
%
%                               Y = X + A^{1/2} * G,                    (1)
%
%   where X is a low-rank signal matrix, G is white noise of variance
%   1/n, and A is an unknown covariance matrix. Pseudo-whitening
%   multiplies by the inverse of the sample covariance
%
%                              \hat A = E*E^T / n0,                     (2)
%
%   where E = A^{1/2} * Z is an independent noise matrix.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors, normalized by sqrt(n)
%   cov0 - the m-by-m sample covariance of noise-only samples
%   m - the dimension
%   n - the number of signal-plus-noise samples
%   gam0 - ratio of dimension to number of noise-only samples
%   k - upper bound on the rank of X
%   iorth - set to 1 to orthogonalize left vectors of \hat{X}
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   xfhat - the estimated matrix, before unwhitening takes place
%   errhat - an estimate of the Frobenius loss between X and \hat{X}
%   khat - estimated rank
%   uhat - approximate left singular vectors of \hat{X} (exact if iorth is 1)
%   vhat - approximate right singular vectors of \hat{X} (exact if iorth is 1)
%   tvals - approximate singular values of \hat{X}
%   etas - generalized singular values of \hat{X}
%   taus - estimated tau parameters
%   sx - estimated singular values of X
%   couts - estimated inner products, on the left
%   cinns - estimated inner products, on the right
%   dmu - estimated average variance
%   
%


%
%        pseudo-whiten the input
%
        [uyfw,uyf,syf,vyf] = fshr_pwhiten(y,cov0,m,n,k);
        dmu = fshr_estmu(cov0,m);

%
%        apply optimal shrinkage
%
        [xhat,xfhat,errhat,khat,uhat,tvals,etas,taus,sx,couts,...
            cinns] = fshr_core(uyfw,uyf,syf,vyf,gam0,dmu,m,n,k,iorth);
        vhat = vyf;

        end
