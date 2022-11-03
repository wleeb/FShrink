        function [xhat0,errhat,topt] = whtd_svdsub(uy,vy,sy,kest,sig,...
            m,n,i1,i2,j1,j2)
%
%                              description:
%
%   Optimal spectral denoiser for submatrix estimation. The input
%   is SVD of the full data matrix, of the form
%
%                               Y = X + G,                          (1)
%
%   where X is an rank k signal matrix, and G is white noise of 
%   variance sig^2. The method attempts to find the optimal k-by-k topt
%   to minimize the AMSE
%
%                || X_0 - \hat U_0 topt \hat V_0^T ||_{Fr}^2        (2)
%
%   where A and B are the diagonal weight matrices, and \hat U_0 and 
%   \hat V_0 are the empirical singular vectors of Y projected onto 
%   coordinates i1 through i2, and j1 through j2, respectively.
%
%   Significantly, this code assumes that the rank estimate kest is correct!
%   No rank estimation is performed by this code.
%
%                               input parameters:
%
%   uy,vy,sy - the svd of the data matrix Y. uy is m-by-k, vy is n-by-k,
%      and sy is k-dimensional vector
%   kest -the rank of the signal matrix X
%   sig - the noise level
%   m,n - the dimensionalities
%   i1,i2 - the first and last rows of X_0
%   j1,j2 - the first and last columnd of X_0
%
%                              output parameters:
%
%   xhat0 - the m2-by-n2 matrix of denoised data
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   topt - estimated k-by-k matrix used in spectral denoiser
%
%
        xhat=zeros(m,n);
        topt=zeros(kest,kest);
        errhat=0;
%
        gam=m/n;
        tol = 100*eps;
%
        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);
%
%        apply weight matrices to empirical vectors
%
        uy2 = uy(i1:i2,:);
        vy2 = vy(j1:j2,:);
%
        amu = (i2-i1+1) / m;
        bmu = (j2-j1+1) / n;
%
%        estimate parameters
%
        dout = uy2'*uy2;
        dinn = vy2'*vy2;
        [tx,eout,einn,coutw,cinnw,cout,cinn] = ...
            whtd_emp2matrs(sy,dout,dinn,kest,amu,bmu,gam,sig);
%
%        solve least squares
%
        topt = whtd_fminprods(tx,kest,dout,dinn,coutw,cinnw,tol);
        [fval,fgrad] = whtd_fevalprods(topt,tx,kest,eout,einn,...
            dout,dinn,coutw,cinnw);
%
%        final estimated submatrix, and estimated error
%
        xhat0 = uy(i1:i2,:) * topt * vy(j1:j2,:)';
        errhat = sqrt(2*fval);

        end
%
