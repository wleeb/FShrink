        function [ell,bedge,err] = mpbdry_thresh(as,bs,awhts,bwhts,m,n,gam)
%
%                            description:
%
%   This code evaluates model parameters for random matrices of the form
%
%                      N = A^{-1/2} G B^{-1/2} / sqrt(n)              (1)
%
%   where G has iid Gaussian entries of variance 1. Returns top eigenvalue
%   of YY^T and the minimum detectable signal eigenvalue. The code is
%   based on the observation that the function F(lambda^*,e) (as function of e)
%   appears to have a double root when lambda^* is the top eigenvalue of YY^T;
%   however, this fact has not been proven in general.
%
%
%                           input parameters:
%
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
%   bedge - estimated eigenvalue of YY^T
%   err - the value of F(bedge,evar); this should be zero
%
%
%
%
%        . . . check input dimensions
%
        if (size(as,1) ~= 1)
%
        as = as';
    end
%
        if (size(bs,1) ~= 1)
%
        bs = bs';
    end
%
        if (size(awhts,1) ~= 1)
%
        awhts = awhts';
    end
%
        if (size(bwhts,1) ~= 1)
%
        bwhts = bwhts';
    end

%
%        . . . evaluate the edge of empirical covariance
%
        bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam);

%
%        evaluate the minimum of F(bedge,e) (as function of e)
%
        evar = mpbdry_fmin(bedge,as,bs,awhts,bwhts,m,n,gam);
        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(evar,bedge,as,bs,awhts,...
            bwhts,m,n,gam);
        eder = -fdlam / fder;
        gderl = gder * eder;
%
%        evaluate the stieltjes transform and derivative
%
        s = sum(1 ./ (gval*as-bedge).*awhts);
        sder = -sum(awhts.*(gderl*as - 1) ./ (gval*as - bedge).^2);
        [sbar,sbar_der] = mpbdry_stra2sbar(s,sder,bedge,gam);
        d = s * sbar * bedge;
%%%        dder = sder*sbar*bedge + s*sbar_der*bedge + s*sbar;

        ell = 1 / d;
%
%        evaluate F(bedge,evar) (should be zero)
%
        ifder=0;
        [err,fder] = mpbdry_evalf(evar,bedge,as,bs,awhts,bwhts,m,n,gam,ifder);


        end
%
