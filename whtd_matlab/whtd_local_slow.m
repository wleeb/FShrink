        function [xhat,errhat] = whtd_local_slow(y,m,n,k,...
            isubr,isubc,nsubr,nsubc,sig)
%
%        performs localized denoising, recomputing SVD for each
%        submatrix
%
%        isubr marks the begining of each segmentation of rows
%        isubc marks the begining of each segmentation of columns
%
        xhat = zeros(m,n);
        isubr = [isubr,m+1];
        isubc = [isubc,n+1];

        errhat = 0;

        for i=1:nsubr
%
        i1 = isubr(i);
        i2 = isubr(i+1) - 1;

        awhts = zeros(1,m);
        awhts(i1:i2) = 1;


        tol = 100*eps;


        for j=1:nsubc
%
        j1 = isubc(j);
        j2 = isubc(j+1) - 1;

        bwhts = zeros(1,n);
        bwhts(j1:j2) = 1;

        [xhat0,errhat0,topt,kest] = whtd_approx(y,sig,awhts,bwhts,...
            m,n,k);
%
        xhat(i1:i2,j1:j2) = xhat0(i1:i2,j1:j2);
        errhat = errhat + errhat0^2;
    end
    end

        errhat = sqrt(errhat);

        end
%
