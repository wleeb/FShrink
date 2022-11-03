        function [ell,cout,cinn] = mpbdry_sback(rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%                            description:
%
%   This code evaluates model parameters for random matrices of the form
%
%                          Y = X + N                               (1)
%
%   where X is a low rank matrix and N = A^{1/2} G B^{1/2} has separable
%   variance profile. Returns eigenvalue of XX^T and angles between
%   singular vectors of X and Y.
%
%                           input parameters:
%
%   rlam - spiked eigenvalue of YY^T
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   ell - estimated eigenvalue of XX^T
%   cout - estimated cosine between left singular vectors of X and Y
%   cinn - estimated cosine between right singular vectors of X and Y
%
%
%        . . . Stieltjes transforms and D transform
%        (derivative of D wrt rlam, NOT sqrt(rlam))
%
        [s,sder] = mpbdry_stiel(rlam,as,bs,awhts,bwhts,m,n,gam);
        [sbar,sbar_der] = mpbdry_stra2sbar(s,sder,rlam,gam);
        d = s * sbar * rlam;
        dder = sder*sbar*rlam + s*sbar_der*rlam + s*sbar;
%
%        spiked model parameters
%
        ell = 1 / d;
        cout = sqrt(s / (dder * ell));
        cinn = sqrt(sbar / (dder * ell));

        end
%
