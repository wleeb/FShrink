        function bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam)
%
%
%                            description:
%
%   This code computes the right boundary of the limiting spectral
%   distribution for a random matrix of the form N = A^{1/2} G B^{1/2},
%   where G is k-by-l with iid entries and the spectrum of A and B
%   are drawn from specified discrete distributions.
%
%   This function's purpose is just to check that the input vectors are
%   correctly oriented (necessary for vectorization) and to call
%   mpbdry_edge0, which does all the actual work.
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
%   bedge - the boundary of the limiting spectral distribution
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

        bedge = mpbdry_edge0(as,bs,awhts,bwhts,m,n,gam);

        end
%
