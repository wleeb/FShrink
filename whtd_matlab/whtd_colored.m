        function [xhat,errhat,topt,kest] = whtd_colored(y,vrows,vcols,...
            m,n,k)
%
%                              description:
%
%   Optimal spectral denoiser for matrix denoising with heteroscedastic
%   noise. The input is a matrix of the form
%
%                               Y = X + G,                           (1)
%
%   where G has noise with rank 1 variance profile and X is a
%   rank k signal matrix. The method whitens the backprojected noise
%   matrix, applies weighted spectral denoising, and applies the inverse
%   transformation.
%
%                               input parameters:
%
%   yholes - backprojected data matrix of size m-by-n
%   vrows,vcols - noise variances of rows and columns; must be COLUMN vectors,
%      of length m and n, respectively
%   m,n - the dimensionalities
%   k - upper bound on the rank of the signal matrix X
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the estimated rank of X
%
%
%
%        . . . whiten the noise
%
        wrows = 1./sqrt(vrows);
        wcols = 1./sqrt(vcols);
        yw = whtd_dmult(y,wrows,wcols,m,n);

%
%        denoise and renormalize
%
        sig=1;
        awhts = sqrt(vrows);
        bwhts = sqrt(vcols);
        [xhat,errhat,topt,kest] = whtd_approx(yw,sig,awhts,bwhts,m,n,k);
        xhat = whtd_dmult(xhat,awhts,bwhts,m,n);

        end
%
