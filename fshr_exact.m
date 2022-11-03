        function [xhat,xfhat,err,tvals,etas] = fshr_exact(y,x,...
            cov0,whts,m,n,n0,k)
%
%                                description:
%
%   Computes exact optimal singular value shrinkage with pseudo-whitening,
%   for input of the form
%
%                               Y = X + A^{1/2} * G,                    (1)
%
%   where X is a low-rank signal matrix, G is white noise of variance
%   1/n, and A is an unknown covariance matrix. Pseudo-whitening
%   multiplies by the inverse of the sample covariance
%
%                              \hat A = E*E^T / n0,                     (2)
%
%   where E = A^{1/2} * Z is an independent noise matrix. Since this code
%   uses the target matrix X, it can be regarded as a ``cheat'' code, useful
%   for theoretical purposes only.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   y - the m-by-n signalmatrix of rank k
%   ep0 - the m-by-n0 noise matrix
%   m - the dimension
%   n - the number of signal-plus-noise samples
%   n0 - the number of noise-only samples
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   xfhat - the estimated matrix, before unwhitening takes place
%   err - the Frobenius loss between X and \hat{X}
%   tvals - approximate singular values of \hat{X}
%   etas - approximate generalized singular values of \hat{X}
%   
%
%

%
        tvals = zeros(k,1);

        cinv0=inv(cov0);
        wmat0=cinv0^(1/2);
        whts=cov0^(1/2);

        yf = wmat0 * y;
        [uyf,syf,vyf] = fshr_svds(yf,m,n,k);
        uyfw = whts*uyf;

%
        [ux,sx,vx] = fshr_svds(x,m,n,k);
%%%        chk0 = norm(ux*diag(sx)*vx' - x,'fro');
%%%        prin2('chk0=',chk0,1);

%%%        prinstop;


%
%        find optimal k-by-k matrix
%
        dout = uyfw'*uyfw;
        dinn = eye(k);
%
        cout = uyfw'*ux;
        cinn = vyf'*vx;

        prinar2('cinn=',cinn,k,k);
        prinar2('cout=',cout,k,k);
        prinar2('dout=',dout,k,k);

%%%        prinstop;
%
        eout = eye(k);
        einn = eye(k);
%
        tol = 100*eps;
        etas = whtd_ddminprods(dout,dinn,cout,cinn,sx,m,n,k,tol);

        prin2('etas=',etas,k);

        [fval,fgrad] = whtd_ddprods(etas,sx,dout,dinn,cout,cinn,...
            eout,einn,m,n,k);

%
%        reconstruct full matrices
%
%%%        xfhat = uyf * diag(etas) * vyf';

        xfhat = fshr_svd2mat(uyf,vyf,etas,m,n,k);
        xhat = whts*xfhat;

        err = sqrt(2*fval);
%%%        tvals = svds(xhat,k);

        for i=1:k
%
        unorm = norm(uyfw(:,i));
        tvals(i) = etas(i) * unorm;
    end


        return;

        [uu,tvals2,vv] = fshr_svds(xhat,m,n,k);

        chk0 = norm(uu*diag(tvals)*vv' - xhat);
        prin2('chk0=',chk0,1);


        xf = wmat0 * x;
        err2 = norm(whts*(xfhat - xf),'fro');
        chk0 = err - err2;

        prin2('chk0=',chk0,1);


        err2 = norm(xhat - x,'fro');

        chk0 = err2-err;

        prin2('err=',err,1);
        prin2('err2=',err2,1);
        prin2('chk0=',chk0,1);

        prinstop;
        end
