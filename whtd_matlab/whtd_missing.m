        function [xhat,errhat,topt,kest] = whtd_missing(yholes,...
            m,n,k,qr,qc,inds,sig)
%
%                              description:
%
%   Optimal spectral denoiser for matrix denoising with missing entries.
%   The input is a matrix of the form
%
%                               Y = X + G,                           (1)
%
%   with certain specified entries replaced with zeros. Here, X is a
%   rank k signal matrix, and G is white noise of variance sig^2. The
%   method whitens the backprojected noise matrix, applies weighted
%   spectral denoising, and applies the inverse transformation.
%
%                               input parameters:
%
%   yholes - backprojected data matrix of size m-by-n
%   m,n - the dimensionalities
%   k - upper bound on the rank of the signal matrix X
%   qr,qc - sampling probabilities on rows and columns; must be COLUMN vectors,
%      of length m and n, respectively
%   inds - 0/1 matrix of size m-by-n marking observed entries
%   sig - the noise level
%   tol - precision for least squares; recommend setting it to ~1d-14
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   topt - estimated kest-by-kest matrix used in spectral denoiser
%   kest - the estimated rank of X
%
%
        yholes = inds.*yholes;
        [yw,awhts,bwhts] = whtd_misswhit(yholes,qr,qc,m,n);
%
%        apply spectral denoising and renormalize
%
        [xhat,errhat,topt,kest] = whtd_approx(yw,sig,awhts,bwhts,m,n,k);
        xhat = whtd_dmult(xhat,awhts,bwhts,m,n);

        end
%
