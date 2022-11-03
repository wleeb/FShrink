        function [rlam,cout,cinn] = mpbdry_sforw(ell,as,bs,awhts,bwhts,m,n,gam)
%
%                            description:
%
%   This code evaluates model parameters for random matrices of the form
%
%                          Y = X + N                               (1)
%
%   where X is a low rank matrix and N = A^{1/2} G B^{1/2} has separable
%   variance profile. Returns asymptotic eigenvalue of YY^T and angles
%   between singular vectors of X and Y.
%
%                           input parameters:
%
%   ell - eigenvalue of XX^T
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   rlam - estimated eigenvalue of YY^T
%   cout - estimated cosine between left singular vectors of X and Y
%   cinn - estimated cosine between right singular vectors of X and Y
%
%
        mmm=1000;

        rlams = zeros(1,mmm);
        vals = zeros(1,mmm);
        ders = zeros(1,mmm);
%
        bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam);
%
%        solve by Newton
%
        zero_mach = mpbdry_machzero();
        tol=sqrt(zero_mach)/10;
        kstop=0;

        rlams(1)=bedge;
        nsteps=mmm;
%
        for i=1:mmm
%
        [val,vder] = mpbdry_dreciproc(rlams(i),as,bs,awhts,bwhts,m,n,gam);
        vals(i)=val - ell;
        ders(i)=vder;
%
        rlams(i+1) = rlams(i) - vals(i) / ders(i);
%
%        check convergence
%
        rr = (rlams(i+1)-rlams(i)) / rlams(i);
        rr=abs(rr);
%%%        prin2('rr=',rr,1);
        if (rr < tol)
%
        kstop=kstop+1;
    end
        if (kstop==2)
%
        nsteps=i;
        break;
    end
    end

        rlam = rlams(nsteps);
        [ell2,cout,cinn] = mpbdry_sback(rlam,as,bs,awhts,bwhts,m,n,gam);
%%%        chk0 = ell2 - ell

        end
%
