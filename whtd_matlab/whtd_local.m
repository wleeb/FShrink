        function [xhat,errhat,kest] = whtd_local(y,m,n,k,...
            isubr,isubc,nsubr,nsubc,sig)
%
%                              description:
%
%   Denoises matrix of the form Y = X + G, where X is low-rank and G
%   is Gaussian noise of specified variance. Applied localized
%   denoising; row and column indices are specified by the user.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%   isubr - integer vector of length nsubr, containing the starting indices
%      of the rows of each submatrix
%   isubc - integer vector of length nsubc, containing the starting indices
%      of the columns of each submatrix
%   nsubr, nsubc - the number of row and column groups, respectively. Total
%      number of submatrices is nsubr*nsubc
%   sig - the standard deviation of the noise
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - the Frobenius loss between X and \hat X
%   kest - estimated rank of X
%
%
        xhat = zeros(m,n);
        isubr = [isubr,m+1];
        isubc = [isubc,n+1];
        errhat = 0;

        gam=m/n;
        tol = 100*eps;

%
%        . . . SVD of data and rank estimation
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%%%        bmargin = 1/n^(2/3);
        bmargin = 0;
        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);
        if (kest == 0)
            return
        end

        for i=1:nsubr
%
        i1 = isubr(i);
        i2 = isubr(i+1) - 1;

        awhts = zeros(1,m);
        awhts(i1:i2) = 1;

        for j=1:nsubc
%
        j1 = isubc(j);
        j2 = isubc(j+1) - 1;

        bwhts = zeros(1,n);
        bwhts(j1:j2) = 1;

        [xhat0,errhat0,topt] = whtd_svdsub(uy,vy,sy,kest,sig,m,n,...
            i1,i2,j1,j2);
%
        xhat(i1:i2,j1:j2) = xhat0;
        errhat = errhat + errhat0^2;
    end
    end

        errhat = sqrt(errhat);

        end
%
