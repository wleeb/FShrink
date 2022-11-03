        function [topt,errhat] = whtd_svdtopt(uy,sy,vy,sig,...
            awhts,bwhts,m,n,k,tol)
%
%                              description:
%
%   Estimates optimal spectral denoiser for data of the form
%
%                               Y = X + G,                          (1)
%
%   where X is a low-rank signal matrix, and G is white noise of 
%   variance sig^2. The method attempts to find the optimal k-by-k topt
%   to minimize the AMSE
%
%                || A (X - \hat U topt \hat V^T) B^T ||_{Fr}^2      (2)
%
%   where A and B are the diagonal weight matrices, and \hat U and \hat V
%   are the empirical singular vectors of Y.
%
%   This code operates on the SVD of Y. Note that no rank estimation
%   is performed; it is assumed that all singular values exceed the 
%   detection threshold.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - the rank of X
%   tol - precision for least squares; recommend setting it to ~1d-14
%
%                              output parameters:
%
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   
%
        topt=zeros(k,k);
        errhat=0;
        
        gam=m/n;
%
%        apply weight matrices to empirical vectors
%
        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);

        amu = sum(awhts.^2) / m;
        bmu = sum(bwhts.^2) / n;
%
%        estimate parameters
%
        dout = uy2'*uy2;
        dinn = vy2'*vy2;
        [tx,eout,einn,coutw,cinnw,cout,cinn] = ...
            whtd_emp2matrs(sy,dout,dinn,k,amu,bmu,gam,sig);
%
%        solve least squares
%
        topt = whtd_fminprods(tx,k,dout,dinn,coutw,cinnw,tol);
%
%        estimate error
%
        [fval,fgrad] = whtd_fevalprods(topt,tx,k,eout,einn,...
            dout,dinn,coutw,cinnw);
        errhat = sqrt(2*fval);


        end
%
