        function [s,sder] = mpbdry_stiel(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%
%                            description:
%
%   This code evaluates the Stieltjes transform and derivative of the
%   limiting spectral distribution for a random matrix of the form 
%   N = A^{1/2} G B^{1/2}, where G is k-by-l with iid entries and the 
%   spectrum of A and B are drawn from specified discrete distributions.
%
%   This function's purpose is just to check that the input vectors are
%   correctly oriented (necessary for vectorization) and to call
%   mpbdry_stiel0, which does all the actual work.
%
%                           input parameters:
%
%   rlam - the value at which we evaluate s and its derivative
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   s - the Stieltjes transform at rlam
%   sder - the derivative of the Stieltjes transform at rlam
%   evar - the value of e(lambda)
%   eder - the value of e'(lambda)
%   gval - the value of G(lambda)
%   gder - derivative of G wrt e
%   gderl - derivative of G wrt lambda 
%
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

        [s,sder,evar,eder,gval,gder,gderl] = mpbdry_stiel0(rlam,...
            as,bs,awhts,bwhts,m,n,gam);

        end
%
