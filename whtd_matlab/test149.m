        function main
        prini(13,1);
%

        test_dense2;

%%%        prinstop;


        test_exact;

        test_local;
        test_missing;
        test_heter;



        prinstop;
        end
%
%
%
%
%
        function test_dense2
%
        randn(1,310);
        rand(1,310);

        m=1500;
        n=3000;

        m=200;
        n=300;

        k=3;
        tx = zeros(k,1);
        sy = zeros(k,1);

        sig=1.4;

        gam=m/n;

        for i=1:k
%
        tx(i) = sig * sqrt(sqrt(gam)) + k - i + 2;
    end

%%%        ua = randn(m,m);
%%%        ua = whtd_gramschmidt(ua,m,m);

        ua=eye(m);

        dd = ua'*ua - eye(m);
        chk0 = norm(dd,'fro');
%%%        prinar2('dd=',dd,m,m);
        prin2('chk0=',chk0,1);

%%%        prinstop;

        awhts = rand(m,1) + .2;
        bwhts = ones(n,1);

        awhts = sort(awhts,'descend');

        adense = ua*diag(awhts)*ua';
        tol=1d-12;


        [x,y,ep,u,v,sx] = genspiked(m,n,k,sig);
        xrot=ua*x;
        yrot=ua*y;

        [xhat2,topt2,errhat2,kest2] = whtd_dense(yrot,sig,adense,m,n,k);

        [xhat3,errhat3,topt3,kest3] = whtd_approx(yrot,sig,awhts,bwhts,m,n,k);
        fff=xhat2-xhat3;
        chk0 = norm(fff,'fro');
        prin2('fff=',fff,10);
        prin2('chk0=',chk0,1);

        end
%
%
%
%
%

        function test_dense
%
        randn(1,310);
        rand(1,310);

        m=1500;
        n=3000;

        m=200;
        n=300;
        k=3;
        tx = zeros(k,1);
        sy = zeros(k,1);

        sig=1.4;

        gam=m/n;

        for i=1:k
%
        tx(i) = sig * sqrt(sqrt(gam)) + k - i + 2;
    end

        ua = randn(m,m);
        ua = whtd_gramschmidt(ua,m,m);

%%%        ua=eye(m);

        dd = ua'*ua - eye(m);
        chk0 = norm(dd,'fro');
%%%        prinar2('dd=',dd,m,m);
        prin2('chk0=',chk0,1);

%%%        prinstop;

        awhts = rand(m,1) + .2;
        bwhts = ones(n,1);

        awhts = sort(awhts,'descend');

        [x,y,ux,vx,uy,vy,sy,amus,bmus,eout,einn,...
            coutw,cinnw,dout,dinn] = whtd_drawcheat(tx,awhts,bwhts,sig,m,n,k);

%%%        prin2('x=',x,m*n);
%%%        prin2('y=',y,m*n);

        tol=1d-12;

        adense = ua*diag(awhts)*ua';

        xrot=ua*x;
        yrot=ua*y;

        k2=k+2;
        [xhat,topt,errhat,kest] = whtd_densecheat(yrot,sig,adense,amus,m,n,k2);
        prin2('xhat=',xhat,m*n);

        dd = adense*(xhat-xrot);

        err=norm(dd,'fro');

        prin2('err=',err,1);
        prin2('errhat=',errhat,1);

        chk0 = err-errhat;
        prin2('chk0=',chk0,1);



        prinstop;
        end
%
%
%
%
%
        function test_exact

        randn(1,200);
        rand(1,100);

        m=40;
        n=20;
        k=3;
        tx = zeros(k,1);
        sy = zeros(k,1);

        sig=1.4;

        gam=m/n;

        for i=1:k
%
        tx(i) = sig * sqrt(sqrt(gam)) + k - i + 2;
    end


        awhts = rand(m,1).^2;
        bwhts = rand(n,1).^2;

        [x,y,ux,vx,uy,vy,sy,amus,bmus,eout,einn,...
            coutw,cinnw,dout,dinn] = whtd_drawcheat(tx,awhts,bwhts,sig,m,n,k);


        prin2('x=',x,m*n);
        prin2('y=',y,m*n);

        tol=1d-12;

        [xhat,uy,vy,errhat,topt,kest] = whtd_cheat(y,sig,...
            awhts,bwhts,amus,bmus,m,n,k);



        [xhat2,errhat2,topt2] = whtd_exact(y,x,awhts,bwhts,m,n,k);

        prin2('topt=',topt,k*k);
        prin2('topt2=',topt2,k*k);


        dif = xhat2 - xhat;
        chk0 = norm(dif,'fro');
        prin2('dif=',dif,12);
        prin2('chk0=',chk0,1);

%%%        prinstop;


%%%        [topt,errhat] = whtd_svdtopt(uy,sy,vy,sig,...
%%%            awhts,bwhts,m,n,k,tol);


        prinar2('topt=',topt,k,k);


        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);
%
        ux2 = whtd_dleft(ux,awhts,m,k);
        vx2 = whtd_dleft(vx,bwhts,n,k);

%%%        ux2(:,2)=-ux2(:,2);
%%%        vx2(:,2)=-vx2(:,2);


        [fval,fgrad] = whtd_fevalvecs(topt,uy2,vy2,ux2,vx2,tx,m,n,k);
        chk0 = norm(fgrad,'fro');
        prin2('fgrad=',fgrad,k*k);
        prin2('topt=',topt,k*k);

        prin2('chk0=',chk0,1);


%%%        prinstop;
        end
%
%
%
%
%
        function test_missing
%
        m=2500;
        n=1700;


%%%        m=500;
%%%        n=400;

        prinf('m=',m,1)
        prinf('n=',n,1)



        k=3;
        sig=1;
        prinf('k=',k,1)

        gam = m/n;

%
%        singular values
%
        sx = sqrt(gam) + 2 + rand(1,k)*20;
        sx = sqrt(sx);

        sx = sort(sx,'descend')';

        [x,y,ep,ux,vx,uy,vy,sy] = whtd_draw(sx,m,n,k,sig);




%
%        generate mask
%
        prows=rand(m,1);
        pcols = rand(n,1);



        pmin=.3;
        pmax = .8;

        for i=1:m
%
        prows(i) = pmin * prows(i) + pmax * (1-prows(i));
    end


        pcols = rand(n,1);

        for j=1:n
%
        pcols(j) = pmin * pcols(j) + pmax * (1-pcols(j));
    end


        pmin2 = min(pcols);
        pmax2 = max(pcols);



        prin2('pmin2 =',pmin2,1)
        prin2('pmax2 =',pmax2,1)

        pmin2 = min(prows);
        pmax2 = max(prows);

        prin2('pmin2 =',pmin2,1)
        prin2('pmax2 =',pmax2,1)

        ps = prows * pcols';


        inds = rand_inds(m,n,ps);


        prinarf('sampled mask=',inds,5,10)

%
%        check sampling probabilities
%
        i=566;
        ppp = sum(inds(i,:)) / n;
        ppp2=prows(i) * sum(pcols) / n;

        prin2('ppp=',ppp,1);
        prin2('ppp2=',ppp2,1);

        j=50;
        ppp=sum(inds(:,j)) / m;
        ppp2=pcols(j) * sum(prows) / m;

        prin2('ppp=',ppp,1);
        prin2('ppp2=',ppp2,1);



        ymask = y .* inds;



        tol=1d-10;

        [xhat,errhat,topt,kest] = whtd_missing(ymask,...
            m,n,k,prows,pcols,inds,sig);

%%%        iii = randi(m*n,1,100);
%%%        xhat(iii)


        err = norm(xhat - x,'fro');

        rel_err = err / norm(sx);

        prin2('err=',err,1)
        prin2('errhat=',errhat,1)

        prin2('rel_err=',rel_err,1)

%%%        prinstop;


        end
%
%
%
%
%
        function test_heter
%
        m=1000;
        n=1500;


%%%        m=500;
%%%        n=400;

        gam = m/n;
        k=3;


        prinf('m=',m,1)
        prinf('n=',n,1)
        prinf('k=',k,1)



%
%        singular values
%
        sx = sqrt(gam) + 1 + rand(1,k)*20;
        sx = sqrt(sx);


        sx(k) = sqrt(gam)-.1;

        sx = sort(sx,'descend');


        prin2('sx=',sx,k)




        mhalf = floor(m/2);
        nhalf = floor(n/2);


%
%        signal singular vectors
%
        u = randn(m,k);

        u(mhalf:m,k) = 0;
        u(1:mhalf,1) = 1;

        u = gramschmidt77(u,m,k);
%
        v(nhalf:n,2) = 0;
        v(1:nhalf,1) = 2;

        v = randn(n,k);
        v = gramschmidt77(v,n,k);


        chk0 = norm(u'*u - eye(k),'fro');
        chk0 = chk0 + norm(v'*v - eye(k),'fro');

        prin2('chk0=',chk0,1)

%
%        noise matrix
%
        avars = rand(m,1) + 1;
        bvars = rand(n,1) + 1;

%%%        avars = ones(m,1);
%%%        bvars = ones(n,1);


        ep = randn(m,n);
        ep = diag(sqrt(avars)) * ep * diag(sqrt(bvars));



%
%        check noise variances correct
%
        i=mhalf+14;
        ddd=sum(ep(i,:).^2) / n;
        ddd2=avars(i)*sum(bvars) / n;

        prin2('ddd=',ddd,1)
        prin2('ddd2=',ddd2,1)

        j=2;
        ddd=sum(ep(:,j).^2) / m;
        ddd2=bvars(j)*sum(avars) / m;

        prin2('ddd=',ddd,1)
        prin2('ddd2=',ddd2,1)
%
%        generate signal and spiked matrices
%
        x = u * diag(sx) * v';
        y = x + ep / sqrt(n);

%

        tol=1d-10;
        [xhat,errhat,topt,kest] = whtd_colored(y,avars,bvars,m,n,k);



%%%        iii = randi(m*n,1,100);
%%%        xhat(iii)

        prinf('kest=',kest,1)

%%%        prinstop;

        dx = xhat - x;
        err = norm(dx,'fro');


        rel_err = err / norm(sx);


        prin2('err=',err,1)
        prin2('errhat=',errhat,1)

        prin2('rel_err=',rel_err,1)

        end
%
%
%
%
%
        function test_local
%
        m=2900;
        n=2600;

        m=500;
        n=400;

        gam = m/n;
        k=4;

        prinf('m=',m,1)
        prinf('n=',n,1)
        prinf('k=',k,1)

        sig=1;

        [x,y,ep,u,v,sx] = genspiked(m,n,k,sig);

%
%        weights are just coordinate projection
%


        m0=floor(m/2);
        n0=floor(n/2) + 10;


        awhts=zeros(m,1);
        awhts(1:m0) = 1;

        bwhts=zeros(n,1);
        bwhts(1:n0) = 1;

        amu = sum(awhts) / m;
        bmu = sum(bwhts) / n;


        isubr = [1,m0+1];
        isubc = [1,n0+1];


        prinf('isubr=',isubr,2)
        prinf('isubc=',isubc,2)



        nsubr = 2;
        nsubc = 2;
        tol=1d-10;
        [xhat,errhat,kest] = whtd_local(y,m,n,k,...
            isubr,isubc,nsubr,nsubc,sig);


%%%        iii = randi(m*n,1,100);
%%%        xhat(iii)

        err = norm(xhat - x,'fro');


        rel_err = err / norm(sx);

        prin2('err=',err,1)
        prin2('errhat=',errhat,1)

        prin2('rel_err=',rel_err,1)

%%%        prinstop

        end
%
%
%
%
%
        function checkpred(tx,coutw,cinnw,dout,dinn,eout,einn,k,amus,bmus,gam,sig)
%
        [sy2,dout2,dinn2,coutw2,cinnw2,cout2,cinn2] = ...
            whtd_pop2cheat(tx,eout,einn,k,amus,bmus,gam,sig);

        dif = dout - dout2;
        chk0=norm(dif,'fro');
        prin2('dif=',dif,k*k);
        prin2('chk0=',chk0,1);
%
        dif = dinn - dinn2;
        chk0=norm(dif,'fro');
        prin2('dif=',dif,k*k);
        prin2('chk0=',chk0,1);
%
        dif = coutw - coutw2;
        chk0=norm(dif,'fro');
        prin2('dif=',dif,k*k);
        prin2('chk0=',chk0,1);
%
        dif = cinnw - cinnw2;
        chk0=norm(dif,'fro');
        prin2('dif=',dif,k*k);
        prin2('chk0=',chk0,1);


        end
%
%
%
%
%
        function checkorth(a,b,whts,m,k)
%
        cc = a'*diag(whts)*b;

        prinar2('cc=',cc,k,k);

%%%        prinstop;


        end
%
%
%
%
%
        function [x,y,ep,u,v,sx] = genspiked(m,n,k,sig)
%
        gam = m/n;

        prin2('gam=',gam,1);
%
%        singular values
%
        sx = sqrt(gam) + 2 + rand(1,k)*20;
        sx = sqrt(sx);

        sx = sort(sx,'descend')';


        for i=1:k
%
        sx(i) = sig * sqrt(sqrt(gam)) + k - i + 2;
    end



        prin2('sx=',sx,k);

%
%        signal singular vectors
%
        u = randn(m,k);
        u = gramschmidt77(u,m,k);
%
        v = randn(n,k);
        v = gramschmidt77(v,n,k);


        chk0 = norm(u'*u - eye(k),'fro');
        chk0 = chk0 + norm(v'*v - eye(k),'fro');

        prin2('chk0=',chk0,1)


        ep = sig*randn(m,n);
%
%        generate signal and spiked matrices
%
        x = u * diag(sx) * v';
        y = x + ep / sqrt(n);



        end
%
%
%
%
%
        function u = gramschmidt77(u,m,k)
%
        u(:,1) = u(:,1) / norm(u(:,1));

        for i=2:k

        for ijk=1:2
%
        for j=1:i-1
%
        pij = sum(u(:,i) .* u(:,j));
        u(:,i) = u(:,i) - pij*u(:,j);
    end
        u(:,i) = u(:,i) / norm(u(:,i));
    end
    end

        end
%
%
%
%
%
        function inds = rand_inds(m,n,ps)
%
        ff = rand(m,n);
        inds = real(ff <= ps);

        end
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       This package contains codes for denoising matrices with weighted
%       loss functions, using spectral denoising. The primary user-callable
%       codes are described below.
%
%   whtd_dense - estimates the optimal spectral denoiser for specified
%       symmetric matrix of row weights.
%
%   whtd_approx - estimates the optimal spectral denoiser for specified
%       row and column weghts.
%
%   whtd_exact - computes the exact spectral denoiser when given the
%       exact low-rank target matrix. For testing purposes only.
%   
%   whtd_local - performs localized denoising, with specified submatrices.
%
%   whtd_missing - performs optimal spectral denoising for matrix wth missing
%      entries. Data is backprojected, noise is whitened, and weighting
%      is performed after denoising.
%
%   whtd_colored - performs optimal spectral denoising for matrix wth
%      additive colored noise with rank 1 variance profile. Noise is
%      whitened, and recoloring is performed after denoising.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
        function [xhat,topt,errhat,kest] = whtd_dense(y,sig,adense,m,n,k)
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
%                || A (X - \hat U topt \hat V^T)  ||_{Fr}^2      (2)
%
%   where A is a symmetric weight matrix, and \hat U and \hat V
%   are the empirical singular vectors of Y.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   adense - m-by-m symmetric weight matrix
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the rank of xhat
%   
%


        xhat=zeros(m,n);
        topt=zeros(k,k);
        errhat=0;
        amus = zeros(k,1);
%
        gam=m/n;
        tol = 100*eps;

%
%        . . . SVD of data
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%        rank estimation
%
        bmargin = 0;
        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);
        if (kest == 0)
%
        return
        end

        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);
%
%        rotate vectors so weights are diagonal
%
        adense = (adense + adense')/2;
        [ua,awhts] = whtd_eigdumb(adense,m);
        uy2=ua'*uy;

%
%        create vectors of weight norms
%
        bwhts = ones(n,1);
        bmus = ones(k,1);
        amu = mean(awhts.^2);
        for i=1:k
%
        amus(i) = amu;
    end
%
%        denoise the matrix
%
        [topt,errhat] = whtd_toptcheat(uy2,sy,vy,sig,...
            awhts,bwhts,amus,bmus,m,n,kest,tol);
%
        xhat = uy * topt * vy';

        end
%
%
%
%
%
        function [xhat,topt,errhat,kest] = whtd_densecheat(y,sig,adense,amus,m,n,k)
%
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
%                || A (X - \hat U topt \hat V^T)  ||_{Fr}^2      (2)
%
%   where A is a symmetric weight matrix, and \hat U and \hat V
%   are the empirical singular vectors of Y.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   adense - m-by-m symmetric weight matrix
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the rank of xhat
%   
%
        xhat=zeros(m,n);
        topt=zeros(k,k);
        errhat=0;

        tol = 100*eps;

%
%        . . . SVD of data
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%        rank estimation
%
        bmargin = 0;
        gam=m/n;
        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);
        if (kest == 0)
%
        return
        end

        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);
%
%        rotate vectors so weights are diagonal
%
        adense = (adense + adense')/2;
        [ua,awhts] = whtd_eigdumb(adense,m);
        uy2=ua'*uy;

%%%        prinf('kest=',kest,1);

        bwhts = ones(n,1);
        bmus = ones(k,1);
        [topt,errhat] = whtd_toptcheat(uy2,sy,vy,sig,...
            awhts,bwhts,amus,bmus,m,n,kest,tol);
%
%        the estimated matrix
%
        xhat = uy * topt * vy';

        end
%
%
%
%
%
        function [xhat,errhat,topt,kest] = whtd_approx(y,sig,...
            awhts,bwhts,m,n,k)
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
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the rank of xhat
%   
%
        xhat=zeros(m,n);
        topt=zeros(k,k);
        errhat=0;

        gam=m/n;
        tol = 100*eps;

%
%        . . . SVD of data
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%        rank estimation
%
        bmargin = 0;
%%%        bmargin = 10/n^(2/3);

        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);
%%%        prinf('kest=',kest,1);


        if (kest == 0)
%
        return
        end
%
        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);

        [topt,errhat] = whtd_svdtopt(uy,sy,vy,sig,...
            awhts,bwhts,m,n,kest,tol);
%
%        the estimated matrix
%
        xhat = uy * topt * vy';

        end
%
%
%
%
%
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
%
%
%
%
        function [xhat,uy,vy,errhat,topt,kest] = whtd_cheat(y,sig,...
            awhts,bwhts,amus,bmus,m,n,k)
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
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the rank of xhat
%   
%
        xhat=zeros(m,n);
        topt=zeros(k,k);
        errhat=0;

        gam=m/n;
        tol = 100*eps;

%
%        . . . SVD of data
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%        rank estimation
%
%%%        bmargin = 1/n^(2/3);
        bmargin = 0;
        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);

        if (kest == 0)
%
        return
        end
%

%%%        prinf('kest=',kest,1);

        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);

        [topt,errhat] = whtd_toptcheat(uy,sy,vy,sig,...
            awhts,bwhts,amus,bmus,m,n,kest,tol);
%
%        the estimated matrix
%
        xhat = uy * topt * vy';

        end
%
%
%
%
%
        function [topt,errhat] = whtd_toptcheat(uy,sy,vy,sig,...
            awhts,bwhts,amus,bmus,m,n,k,tol)
%
        topt=zeros(k,k);
        errhat=0;
        
        gam=m/n;
%
%        apply weight matrices to empirical vectors
%
        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);
%
%        estimate parameters
%
        dout = uy2'*uy2;
        dinn = vy2'*vy2;
        [tx,eout,einn,coutw,cinnw,cout,cinn] = ...
            whtd_emp2cheat(sy,dout,dinn,k,amus,bmus,gam,sig);

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
%
%
%
%
        function [ts,eout,einn,coutw,cinnw,cout,cinn] = ...
            whtd_emp2cheat(sy,dout,dinn,k,amus,bmus,gam,sig)
%
        ells=zeros(k,1);
        cout=zeros(k,1);
        cinn=zeros(k,1);
        ts=zeros(k,1);
%
%        parameters for each individual component
%
        for i=1:k
%
        [ts(i),cout(i),cinn(i)] = whtd_emp2pop(sy(i),gam,sig);
    end
%
%        build E matrices of weighted inner products
%
        eout = diag(1./cout) * dout * diag(1./cout);
        einn = diag(1./cinn) * dinn * diag(1./cinn);
        for i=1:k
%
        eout(i,i) = (dout(i,i) - (1-cout(i)^2)*amus(i)) / cout(i)^2;
        einn(i,i) = (dinn(i,i) - (1-cinn(i)^2)*bmus(i)) / cinn(i)^2;
    end
%
%        build C matrices of weighted inner products
%
        coutw = diag(cout) * eout;
        cinnw = diag(cinn) * einn;

        end
%
%
%
%
%
        function [xhat,errhat,topt] = whtd_exact(y,x,awhts,bwhts,m,n,k)
%
%                              description:
%
%   Evaluates optimal k-by-k matrix topt for minimizing
%
%                || A (X - \hat U topt \hat V^T) B^T ||_{Fr}^2      (1)
%
%   where A and B are the diagonal weight matrices, and \hat U and \hat V
%   are the empirical singular vectors of Y.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   x - the m-by-n matrix of data vectors
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - the Frobenius loss between X and \hat X
%   topt - optimal k-by-k matrix used in spectral denoiser
%   
%
%        . . . apply weight matrices to empirical vectors
%

        tol = 100*eps;
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);

%
%        apply weight matrices to population vectors
%
        [ux,sx,vx] = whtd_svdsmart(x,m,n,k);
        ux2 = whtd_dleft(ux,awhts,m,k);
        vx2 = whtd_dleft(vx,bwhts,n,k);

%
%        solve least squares for k-by-k topt
%
        topt = whtd_fminvecs(uy2,vy2,ux2,vx2,sx,m,n,k,tol);
        [fval,fgrad] = whtd_fevalvecs(topt,uy2,vy2,ux2,vx2,sx,m,n,k);
%%%        chk0 = norm(uy2 * topt * vy2' - ux2 * diag(sx) * vx2','fro') ...
%%%            - sqrt(2*fval)
%
%        final matrix and error
%
        xhat = uy * topt * vy';
        errhat = sqrt(2*fval);
%%%        err = norm(diag(awhts)*(xhat - x)*diag(bwhts),'fro');
%%%        chk0 = errhat - err

        end
%
%
%
%
%
        function [xhat,errhat,topt,kest] = whtd_ddapprox(y,sig,...
            awhts,bwhts,m,n,k)
%
%                              description:
%
%   Estimates optimal spectral denoiser for data of the form
%
%                               Y = X + G,                          (1)
%
%   where X is a low-rank signal matrix, and G is white noise of 
%   variance sig^2. The method attempts to find the optimal k-by-k 
%   diagonal matrix topt to minimize the AMSE
%
%                || A (X - \hat U topt \hat V^T) B^T ||_{Fr}^2      (2)
%
%   where A and B are the diagonal weight matrices, and \hat U and \hat V
%   are the empirical singular vectors of Y. This code differs from
%   whtd_approx since topt is required to be diagonal.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the rank of xhat
%   
%
        xhat=zeros(m,n);
        topt=zeros(k,k);
        errhat=0;
        kest=0;

        gam=m/n;
        tol = 100*eps;

%
%        . . . SVD of data
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%        rank estimation
%
%%%        bmargin = 1/n^(2/3);
        bmargin = 0;
        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);

        if (kest == 0)
%
        return
        end
%
        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);

        [topt,errhat] = whtd_ddtopt(uy,sy,vy,sig,...
            awhts,bwhts,m,n,kest,tol);
%
%        the estimated matrix
%
        xhat = uy * diag(topt) * vy';

        end
%
%
%
%
%
        function [topt,errhat] = whtd_ddtopt(uy,sy,vy,sig,...
            awhts,bwhts,m,n,k,tol)
%
        topt=zeros(k,1);
        errhat=0;

        amu = sum(awhts.^2) / m;
        bmu = sum(bwhts.^2) / n;

        prin2('amu=',amu,1);

        gam=m/n;
%
%        apply weight matrices to empirical vectors
%
        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);
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
        topt = whtd_ddminprods(dout,dinn,coutw,cinnw,tx,m,n,k,tol);
%
%        estimate error
%
        [fval,fgrad] = whtd_ddprods(topt,tx,dout,dinn,coutw,cinnw,...
            eout,einn,m,n,k);
        errhat = sqrt(2*fval);



        end
%
%
%
%
%
        function [xhat,errhat,topt,kest] = whtd_ddcheat(y,sig,...
            awhts,bwhts,amus,bmus,m,n,k)
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
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the rank of xhat
%   
%
        xhat=zeros(m,n);
        topt=zeros(k,k);
        errhat=0;

        gam=m/n;
        tol = 100*eps;

%
%        SVD of data
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%        rank estimation
%
%%%        bmargin = 1/n^(2/3);
        bmargin = 0;
        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);

        if (kest == 0)
%
        return
        end
%
        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);
%
%        optimal coefficients
%
        [topt,errhat] = whtd_ddtoptcheat(uy,sy,vy,sig,...
            awhts,bwhts,amus,bmus,m,n,kest,tol);
%
%        the estimated matrix
%
        xhat = uy * diag(topt) * vy';

        end
%
%
%
%
%
        function [topt,errhat] = whtd_ddtoptcheat(uy,sy,vy,sig,...
            awhts,bwhts,amus,bmus,m,n,k,tol)
%
        topt=zeros(k,1);
        errhat=0;
        
        gam=m/n;
%
%        apply weight matrices to empirical vectors
%
        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);
%
%        estimate parameters
%
        dout = uy2'*uy2;
        dinn = vy2'*vy2;
        [tx,eout,einn,coutw,cinnw,cout,cinn] = ...
            whtd_emp2cheat(sy,dout,dinn,k,amus,bmus,gam,sig);
%
%        solve least squares
%
        topt = whtd_ddminprods(dout,dinn,coutw,cinnw,tx,m,n,k,tol);

%
%        estimate error
%
        [fval,fgrad] = whtd_ddprods(topt,tx,dout,dinn,coutw,cinnw,...
            eout,einn,m,n,k);
        errhat = sqrt(2*fval);



        end
%
%
%
%
%
        function [xhat,errhat,topt] = whtd_ddexact(y,x,awhts,bwhts,m,n,k)
%
%                              description:
%
%   Evaluates optimal k-by-k diagonal matrix topt for minimizing
%
%                || A (X - \hat U topt \hat V^T) B^T ||_{Fr}^2      (1)
%
%   where A and B are the diagonal weight matrices, and \hat U and \hat V
%   are the empirical singular vectors of Y.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   x - the m-by-n matrix of data vectors
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - the Frobenius loss between X and \hat X
%   topt - optimal k-by-k matrix used in spectral denoiser
%   
%
%        . . . apply weight matrices to empirical vectors
%

        tol = 100*eps;
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
        uy2 = whtd_dleft(uy,awhts,m,k);
        vy2 = whtd_dleft(vy,bwhts,n,k);
%
%        apply weight matrices to population vectors
%
        [ux,tx,vx] = whtd_svdsmart(x,m,n,k);
        ux2 = whtd_dleft(ux,awhts,m,k);
        vx2 = whtd_dleft(vx,bwhts,n,k);
%
%        solve least squares for k-by-k topt
%
        topt = whtd_ddminvecs(uy2,vy2,tx,ux2,vx2,m,n,k,tol);
        [fval,fgrad] = whtd_ddvecs(topt,uy2,vy2,tx,ux2,vx2,m,n,k);
%
%        final matrix and error
%
        xhat = uy * diag(topt) * vy';
        errhat = sqrt(2*fval);

        end
%
%
%
%
%
        function tmin = whtd_ddminprods(dout,dinn,cout,cinn,tx,m,n,k,tol)
%
        aa = dout.*dinn;
        bb = cout.*cinn;

        [ainv,khat] = whtd_pseudoinv(aa,k,k,tol);
        tmin = ainv * bb * tx;

        end
%
%
%
%
%
        function [fval,fgrad] = whtd_ddprods(tvec,tx,dout,dinn,cout,cinn,...
            eout,einn,m,n,k)
%
        fgrad = zeros(k,1);

        prin2('tvec=',tvec,k);
        prin2('tx=',tx,k);

        aa1=dout * diag(tvec) * dinn;
        prod1 = whtd_scapro(diag(tvec),aa1,k*k);
%
        aa2=eout * diag(tx) * einn;
        prod2 = whtd_scapro(diag(tx),aa2,k*k);
%
        aa3=cout*diag(tx)*cinn';
        prod3 = whtd_scapro(diag(tvec),aa3,k*k);
%
        fval = prod1 + prod2 - 2*prod3;
        fval = fval/2;


        aa = dout.*dinn;
        bb = cout.*cinn;
        fgrad = aa*tvec - bb*tx;

        end
%
%
%
%
%
        function tmin = whtd_ddminvecs(uy,vy,tx,ux,vx,m,n,k,tol)
%
        dout = uy'*uy;
        dinn = vy'*vy;
        coutw = uy'*ux;
        cinnw = vy'*vx;
%
        aa = dout.*dinn;
        bb = coutw.*cinnw;

        [ainv,khat] = whtd_pseudoinv(aa,k,k,tol);
        tmin = ainv * bb * tx;

        end
%
%
%
%
%
        function [fval,fgrad] = whtd_ddvecs(tvec,uy,vy,tx,ux,vx,m,n,k)
%
        prin2('tvec=',tvec,k);
        prin2('tx=',tx,k);
%
        ahat = uy*diag(tvec)*vy';
        arhs = ux*diag(tx)*vx';
        adif = ahat - arhs;
        fval = norm(adif,'fro')^2 / 2;

%
%        evaluate the gradient of f
%
        fgrad = zeros(k,1);

        dout = uy'*uy;
        dinn = vy'*vy;
        coutw = uy'*ux;
        cinnw = vy'*vx;
%
        aa = dout.*dinn;
        bb = coutw.*cinnw;
        fgrad = aa*tvec - bb*tx;

        end
%
%
%
%
%
        function [xhat,errhat,kest] = whtd_local(y,m,n,k,...
            isubr,isubc,nsubr,nsubc,sig)
%
%                              description:
%
%   Denoises matrix of the form Y = X + G, where X is low-rank and G
%   is Gaussian noise of specified variance. Applied localized
%   denoising; row and column indices are specified by the user.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%   isubr - integer vector of length nsubr, containing the starting indices
%      of the rows of each submatrix
%   isubc - integer vector of length nsubc, containing the starting indices
%      of the columns of each submatrix
%   nsubr, nsubc - the number of row and column groups, respectively. Total
%      number of submatrices is nsubr*nsubc
%   sig - the standard deviation of the noise
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - the Frobenius loss between X and \hat X
%   kest - estimated rank of X
%
%
        xhat = zeros(m,n);
        isubr = [isubr,m+1];
        isubc = [isubc,n+1];
        errhat = 0;

        gam=m/n;
        tol = 100*eps;

%
%        . . . SVD of data and rank estimation
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%%%        bmargin = 1/n^(2/3);
        bmargin = 0;
        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);
        if (kest == 0)
            return
        end

        for i=1:nsubr
%
        i1 = isubr(i);
        i2 = isubr(i+1) - 1;

        awhts = zeros(1,m);
        awhts(i1:i2) = 1;

        for j=1:nsubc
%
        j1 = isubc(j);
        j2 = isubc(j+1) - 1;

        bwhts = zeros(1,n);
        bwhts(j1:j2) = 1;

        [xhat0,errhat0,topt] = whtd_svdsub(uy,vy,sy,kest,sig,m,n,...
            i1,i2,j1,j2);
%
        xhat(i1:i2,j1:j2) = xhat0;
        errhat = errhat + errhat0^2;
    end
    end

        errhat = sqrt(errhat);

        end
%
%
%
%
%
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
%
%
%
%
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
%
%
%
%
        function [yw,awhts,bwhts] = whtd_misswhit(yholes,qr,qc,m,n)
%
        awhts = sqrt(1./qr);
        bwhts = sqrt(1./qc);
%
%        equalize the noise variances
%
        yw = whtd_dmult(yholes,awhts,bwhts,m,n);

        end
%
%
%
%
%
        function a2 = whtd_dmult(a,dl,dr,m,n)
%
%        forms the product
%
%             A2 = Dleft * A * Dright,                        (1)
%
%        where Dleft and Dright are specified diagonal matrices
%
%
%        . . . ensure dl and dr are column vectors
%
        [n1,n2] = size(dl);
        if (n1 < n2)
%
        dl=dl';
    end
        [n1,n2] = size(dr);
        if (n1 < n2)
%
        dr=dr';
    end
%
%        form the product
%
        a2 = bsxfun(@times,a,dl);
        a2 = bsxfun(@times,a2',dr)';

        end
%
%
%
%
%
        function a2 = whtd_dleft(a,d,m,n)
%
%        forms the product
%
%             A2 = D * A,                               (1)
%
%        where D is a specified diagonal matrix
%

%
%        ensure diag is a column vector
%
        [n1,n2] = size(d);
        if (n1 < n2)
%
        d=d';
    end
%
%        form the product
%
        a2 = bsxfun(@times,a,d);

        end
%
%
%
%
%
        function [xhat,errhat,topt,kest] = whtd_colored(y,vrows,vcols,...
            m,n,k)
%
%                              description:
%
%   Optimal spectral denoiser for matrix denoising with heteroscedastic
%   noise. The input is a matrix of the form
%
%                               Y = X + G,                           (1)
%
%   where G has noise with rank 1 variance profile and X is a
%   rank k signal matrix. The method whitens the backprojected noise
%   matrix, applies weighted spectral denoising, and applies the inverse
%   transformation.
%
%                               input parameters:
%
%   yholes - backprojected data matrix of size m-by-n
%   vrows,vcols - noise variances of rows and columns; must be COLUMN vectors,
%      of length m and n, respectively
%   m,n - the dimensionalities
%   k - upper bound on the rank of the signal matrix X
%
%                              output parameters:
%
%   xhat - the m-by-n matrix of denoised data
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the estimated rank of X
%
%
%
%        . . . whiten the noise
%
        wrows = 1./sqrt(vrows);
        wcols = 1./sqrt(vcols);
        yw = whtd_dmult(y,wrows,wcols,m,n);

%
%        denoise and renormalize
%
        sig=1;
        awhts = sqrt(vrows);
        bwhts = sqrt(vcols);
        [xhat,errhat,topt,kest] = whtd_approx(yw,sig,awhts,bwhts,m,n,k);
        xhat = whtd_dmult(xhat,awhts,bwhts,m,n);

        end
%
%
%
%
%
        function [ts,eout,einn,coutw,cinnw,cout,cinn] = ...
            whtd_emp2matrs(sy,dout,dinn,k,amu,bmu,gam,sig)
%
        ells=zeros(k,1);
        cout=zeros(k,1);
        cinn=zeros(k,1);
        ts=zeros(k,1);
%
%        parameters for each individual component
%
        for i=1:k
%
        [ts(i),cout(i),cinn(i)] = whtd_emp2pop(sy(i),gam,sig);
    end
%
%        build E matrices of weighted inner products
%
        eout = diag(1./cout) * dout * diag(1./cout);
        einn = diag(1./cinn) * dinn * diag(1./cinn);
        for i=1:k
%
        eout(i,i) = (dout(i,i) - (1-cout(i)^2)*amu) / cout(i)^2;
        einn(i,i) = (dinn(i,i) - (1-cinn(i)^2)*bmu) / cinn(i)^2;
    end
%
%        build C matrices of weighted inner products
%
        coutw = diag(cout) * eout;
        cinnw = diag(cinn) * einn;

        end
%
%
%
%
%
        function [tx,cout,cinn] = whtd_emp2pop(rlam,gam,sig)
%
        rlam=rlam/sig;

        btop = 1+sqrt(gam);
        if (rlam <= btop);
%
        tx=0;
        cout=0;
        cinn=0;
        return;
    end

        yy = rlam^2 - 1 - gam;
        tx = (yy + sqrt(yy^2 - 4*gam))/2;
        tx = sqrt(tx);
%
        cout = (tx^4 - gam) / (tx^4 + tx^2*gam);
        cout = sqrt(cout);
%
        cinn = (tx^4 - gam) / (tx^4 + tx^2);
        cinn = sqrt(cinn);

        tx=tx*sig;

        end
%
%
%
%
%
        function [sy,dout,dinn,coutw,cinnw,cout,cinn] = ...
            whtd_pop2matrs(tx,eout,einn,k,amu,bmu,gam,sig);
%
        sy=zeros(1,k);
        cout=zeros(k,1);
        cinn=zeros(k,1);
        coutw=zeros(k,k);
        cinnw=zeros(k,k);
        dout=zeros(k,k);
        dinn=zeros(k,k);
%
%        parameters for each individual component
%
        for i=1:k
%
        [sy(i),cout(i),cinn(i)] = whtd_pop2emp(tx(i),gam,sig);
    end
%
%        build D matrices of weighted inner products
%
        dout = diag(cout) * eout * diag(cout);
        dinn = diag(cinn) * einn * diag(cinn);
        for i=1:k
%
        dout(i,i) = cout(i)^2*eout(i,i) +  (1-cout(i)^2)*amu;
        dinn(i,i) = cinn(i)^2*einn(i,i) +  (1-cinn(i)^2)*bmu;
    end
%
%        build C matrices of weighted inner products
%
        coutw = diag(cout) * eout;
        cinnw = diag(cinn) * einn;

        end
%
%
%
%
%
        function [sy,cout,cinn] = whtd_pop2emp(tx,gam,sig)
%
        tx=tx/sig;
%
        tmin = sqrt(sqrt(gam));
        if (tx <= tmin);
%
        sy=sig*(1 + sqrt(gam));
        cout=0;
        cinn=0;
        return;
    end

        sy = (tx^2 + 1)*(1 + gam/tx^2);
        sy = sig*sqrt(sy);

        cout = (tx^4 - gam) / (tx^4 + tx^2*gam);
        cout = sqrt(cout);

        cinn = (tx^4 - gam) / (tx^4 + tx^2);
        cinn = sqrt(cinn);

        end
%
%
%
%
%
        function [sy,dout,dinn,coutw,cinnw,cout,cinn] = ...
            whtd_pop2cheat(tx,eout,einn,k,amus,bmus,gam,sig);
%
        sy=zeros(k,1);
        cout=zeros(k,1);
        cinn=zeros(k,1);
        coutw=zeros(k,k);
        cinnw=zeros(k,k);
        dout=zeros(k,k);
        dinn=zeros(k,k);
%
%        parameters for each individual component
%
        for i=1:k
%
        [sy(i),cout(i),cinn(i)] = whtd_pop2emp(tx(i),gam,sig);
    end
%
%        build D matrices of weighted inner products
%
        dout = diag(cout) * eout * diag(cout);
        dinn = diag(cinn) * einn * diag(cinn);
        for i=1:k
%
        dout(i,i) = cout(i)^2*eout(i,i) +  (1-cout(i)^2)*amus(i);
        dinn(i,i) = cinn(i)^2*einn(i,i) +  (1-cinn(i)^2)*bmus(i);
    end
%
%        build C matrices of weighted inner products
%
        coutw = diag(cout) * eout;
        cinnw = diag(cinn) * einn;

        end
%
%
%
%
%
        function kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin)
%
        rlams = sy.^2;
%
        bedge = sig^2*(1 + sqrt(gam))^2;
        bedge = bedge + sig^2*bmargin;
        kest=0;
        for i=1:k
%
        if (rlams(i) <= bedge)
            break;
        end
        kest = kest+1;
    end

        end
%
%
%
%
%
        function tmin = whtd_fminprods(t,k,dout,dinn,cout,cinn,tol)
%
        [dout_inv,khat] = whtd_pseudoinv(dout,k,k,tol);
        [dinn_inv,khat] = whtd_pseudoinv(dinn,k,k,tol);
        tmin = dout_inv*cout*diag(t)*cinn'*dinn_inv';

        end
%
%
%
%
%
        function [fval,fgrad] = whtd_fevalprods(t_hat,t,k,...
            eout,einn,dout,dinn,cout,cinn)
%
        aa1=dout * t_hat * dinn;
        prod1 = whtd_scapro(t_hat,aa1,k*k);
%
        aa2=eout * diag(t) * einn;
        prod2 = whtd_scapro(diag(t),aa2,k*k);
%
        aa3=cout*diag(t)*cinn';
        prod3 = whtd_scapro(t_hat,aa3,k*k);
%
        fval = prod1 + prod2 - 2*prod3;
        fval = fval/2;

        fgrad = dout * t_hat * dinn - cout*diag(t)*cinn';

        end
%
%
%
%
%
        function prod = whtd_scapro(x,y,n)
%
        prod=sum(x(1:n).*y(1:n));

        return;

        x=x(:);
        y=y(:);

        prod = 0;
        for i=1:n
%
        prod = prod + x(i)*y(i);
    end
        end
%
%
%
%
%
        function tmin = whtd_fminvecs(what,zhat,w,z,t,m,n,k,tol)
%
        [winv,khat] = whtd_pseudoinv(what,m,k,tol);
        [zinv,khat] = whtd_pseudoinv(zhat,n,k,tol);
        tmin = winv*w*diag(t)*z'*zinv';

        end
%
%
%
%
%
        function [fval,fgrad] = whtd_fevalvecs(tvar,what,zhat,w,z,t,m,n,k)
%
        xhat = what * tvar * zhat';
        x = w*diag(t)*z';
%
        fval = norm(xhat - x,'fro')^2 / 2;
        fgrad = what' * (xhat - x) * zhat;

        end
%
%
%
%
%
        function [winv,khat] = whtd_pseudoinv(w,m,n,tol)
%
        khat=min(m,n);
        [wl,ws,wr] = whtd_svdsmart(w,m,n,khat);
%%%        chk0 = norm(wl * diag(ws) * wr' - w,'fro')

        for i=1:min(m,n)
%
        rri = ws(i) / ws(1);
        if (rri < tol)
%
        khat = i - 1;
        break;
    end
    end

        wr = wr(:,1:khat);
        wl = wl(:,1:khat);
        ws = ws(1:khat);
        winv = wr * diag(1./ws) * wl';

        end
%
%
%
%
%
        function [u,s,v] = whtd_svdsmart(a,m,n,k)
%
        if (m/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
    end

        if (m/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
    end

        end
%
%
%
%
%
        function [ua,sa] = whtd_eigdumb(a,m)
%
        [ua,sa]=eig(a);
        sa=diag(sa);
        [sa,ii] = sort(sa,'descend');
        ua=ua(:,ii);

        end
%
%
%
%
%
        function wnorm = whtd_norm(x,awhts,bwhts,m,n)
%
        x2 = whtd_dmult(x,awhts,bwhts,m,n);
        wnorm = norm(x2,'fro');

        end
%
%
%
%
%
        function whtd_mpcompare(y,m,n,sig,ifig)
%
        s = svd(y);
        evals = s.^2;
        gam = m/n;

        whtd_mpcompare0(evals,gam,sig,ifig);

        end
%
%
%
%
%
        function whtd_mpcompare0(evals,gam,sig,ifig)
%
%        plots histogram of evals and overlays plot of continuous
%        Marchenko-Pastur density, white noise with variance 1,
%        with aspect ratio gam. ifig is the figure number of plot
%
        x0 = sig^2*(1-sqrt(gam))^2;
        x1 = sig^2*(1+sqrt(gam))^2;

        nvals=100;
        ts=linspace(x0,x1,nvals);

        for i=1:nvals
%
        vals_mp(i) = whtd_mpeval(ts(i),gam,sig);
    end

        figure(ifig);
        subplot(1,2,1);
        hh = histogram(evals,100,'Normalization','pdf');
        hold on; plot(ts,vals_mp,'linewidth',3)

        m=length(evals);

        subplot(1,2,2);
        plot(evals,'*','Color','b')
        hold on; plot(0:m+1,x1*ones(m+2,1),'LineWidth',2,'Color','r')
        hold on; plot(0:m+1,x0*ones(m+2,1),'LineWidth',2,'Color','r')

        xlim([0,m+1])

        set(figure(ifig),'Position',[500,500,1300,500])

        end
%
%
%
%
%
        function val = whtd_mpeval(t,gam,sig)
%
%        evaluates standard Marchenko-Pastur density with 
%        aspect ratio gam at t
%
        x0 = sig^2*(1-sqrt(gam))^2;
        x1 = sig^2*(1+sqrt(gam))^2;
        val = (x1 - t) * (t - x0);
        val = sqrt(val) / (2*pi*t*gam) / sig^2;

        end
%
%
%
%
%
        function [err,rlam,cout,cinn] = whtd_errfmla(sx,alpha,...
            beta,gam,sig,amu,bmu)
%
%        evaluates weighted Frobenius norm error for rank 1 signal
%
        [rlam,cout,cinn] = whtd_pop2emp(sx,gam,sig);

        sout = sqrt(1 - cout^2);
        sinn = sqrt(1 - cinn^2);

        eta1 = alpha / (cout^2*alpha + sout^2*amu);
        eta2 = beta / (cinn^2*beta + sinn^2*bmu);
        eta = eta1 * eta2;

        err = sx^2*alpha*beta*(1 - cout^2*cinn^2 * eta);
        err = sqrt(err);

        end
%
%
%
%
%
        function [phi,phi_der,rphi] = whtd_evalphi(s,gam,sig,amu,alpha)
%
%        rescale alpha and amu
%
        alpha = alpha / amu;
        amu=1;

        [cout,cout_der,rout] = whtd_evalcout(s,gam,sig);
%
%        phi and derivative phi'
%
        phi = alpha * cout / (alpha*cout^2 + amu*(1-cout^2));
        phi_der = alpha*cout_der*(1 - (alpha-1)*cout^2) ...
            / (1 + (alpha-1)*cout^2)^2;
%
%        ratio between phi' and phi
%
        dd = (alpha-1)*cout^2;
        ff = (1 - dd) / (1 + dd);
        rphi = (cout_der / cout) * ff;

        end
%
%
%
%
%
        function [psi,psi_der,rpsi] = whtd_evalpsi(s,gam,sig,bmu,beta)
%
%        rescale beta and bmu
%
        beta = beta / bmu;
        bmu = 1;

        [cinn,cinn_der,rinn] = whtd_evalcinn(s,gam,sig);

%
%        psi and derivative psi'
%
        psi = beta * cinn / (beta*cinn^2 + bmu*(1-cinn^2));
        psi_der = beta*cinn_der*(1 - (beta-1)*cinn^2) ...
            / (1 + (beta-1)*cinn^2)^2;
%
%        ratio between psi' and psi
%
        dd = (beta-1)*cinn^2;
        ff = (1 - dd) / (1 + dd);
        rpsi = (cinn_der / cinn) * ff;

        end
%
%
%
%
%
        function [topt,topt_der,rtopt] = whtd_evaltopt(s,gam,sig,amu,bmu,...
           alpha,beta)
%
%        evaluate optimal singular value for rank 1 signal, and its
%        derivative, and ratio of derivative to value
%
        [rlam,cout,cinn] = whtd_pop2emp(s,gam,sig);

        [phi,phi_der,rphi] = whtd_evalphi(s,gam,sig,amu,alpha);
        [psi,psi_der,rpsi] = whtd_evalpsi(s,gam,sig,bmu,beta);

        topt = s*phi*psi;
        topt_der = phi*psi + s*phi_der*psi + s*phi*psi_der;
%
%        formula for ratio of derivative to topt
%
        rtopt = phi_der / phi + psi_der / psi + 1/s;

        end
%
%
%
%
%
        function [cout,cout_der,rout] = whtd_evalcout(s,gam,sig)
%
%        evaluates outer cosine cout, its derivative, and the
%        ratio between derivative and the value
%
        s0 = s/sig;
        ell0=s0^2;
        sig0=1;
        [rlam,cout,cinn] = whtd_pop2emp(s0,gam,sig0);
%
        csq = cout^2;
        csq_der = 2*gam*s0*(s0^4 + 2*s0^2 + gam) / (s0^4 + gam*s0^2)^2;
        cout_der = csq_der / cout / 2;
        cout_der = cout_der / sig;
%
%        ratio between derivative and value
%
        rout = gam*(s0^4 + 2*s0^2 + gam) / (s0 * (s0^2+gam) * (s0^4-gam));
        rout = rout / sig;

        end
%
%
%
%
%
        function [cinn,cinn_der,rinn] = whtd_evalcinn(s,gam,sig)
%
%        evaluates inner cosine cinn, its derivative, and the
%        ratio between derivative and the value
%
        s0 = s/sig;
        ell0=s0^2;
        sig0=1;
        [rlam,cout,cinn] = whtd_pop2emp(s0,gam,sig0);
%
        csq = cinn^2;
        csq_der = 2*s0*(s0^4 + 2*gam*s0^2 + gam) / (s0^4 + s0^2)^2;
        cinn_der = csq_der / cinn / 2;

        cinn_der = cinn_der / sig;
%
%        ratio between derivative and value
%
        rinn = (s0^4 + 2*gam*s0^2 + gam) / (s0*(s0^2+1)*(s0^4-gam));
        rinn = rinn / sig;

        end
%
%
%
%
%
        function [x,y,ep,ux,vx,uy,vy,sy] = whtd_draw(sx,m,n,k,sig)
%
        gam = m/n;

%%%        prin2('gam=',gam,1);
%%%        prin2('sx=',sx,k);
%
%        signal singular vectors
%
        ux = randn(m,k);
        ux = whtd_gramschmidt(ux,m,k);
%
        vx = randn(n,k);
        vx = whtd_gramschmidt(vx,n,k);

%%%        checkorth(u,u,ones(m,1),m,k);
%%%        checkorth(v,v,ones(n,1),n,k);
%%%        prinstop;

%
%        add iid Gaussian noise
%
        ep = randn(m,n);
        ep = sig * ep / sqrt(n);
%
%        generate signal and spiked matrices
%
        x = ux * diag(sx) * vx';
        y = x + ep;

        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);

        end
%
%
%
%
%
        function [x,y,ux,vx,uy,vy,sy,amus,bmus,eout,einn,...
            coutw,cinnw,dout,dinn] = whtd_drawcheat(tx,awhts,bwhts,sig,m,n,k)
%
        couts=zeros(k,1);
        cinns=zeros(k,1);
        sy=zeros(k,1);

        awhts = awhts.^2;
        bwhts = bwhts.^2;

        gam=m/n;

%
%        evaluate cosines
%
        for i=1:k
%
        [sy(i),couts(i),cinns(i)] = whtd_pop2emp(tx(i),gam,sig);
    end
%
%        generate vectors
%
        [ux,uy,ures,amus] = whtd_vecs(couts,m,k,awhts,gam,sig);
        [vx,vy,vres,bmus] = whtd_vecs(cinns,n,k,bwhts,gam,sig);

        y = uy * diag(sy) * vy';
        x = ux * diag(tx) * vx';
%
%        inner product matrices
%
        eout = ux'*diag(awhts)*ux;
        einn = vx'*diag(bwhts)*vx;

        dout = uy'*diag(awhts)*uy;
        dinn = vy'*diag(bwhts)*vy;

        coutw = uy'*diag(awhts)*ux;
        cinnw = vy'*diag(bwhts)*vx;

%%%        checkpred(tx,coutw,cinnw,dout,dinn,eout,einn,k,amus,bmus,gam,sig);

        end
%
%
%
%
%
        function [ux,uy,ures,amus] = whtd_vecs(cs,m,k,whts,gam,sig)
%
        ss = zeros(k,1);
        amus = zeros(k,1);
        uy = zeros(m,k);
%
%%%        prinf('m=',m,1);
%%%        prinf('k=',k,1);

%
%        random orthonormal X vectors
%
        ux=randn(m,k);
        ux = whtd_gramschmidt(ux,m,k);
%
%        residual vectors
%
        ures = whtd_res(ux,whts,m,k);
%
%        contruct vectors for Y
%
        for i=1:k
%
        amus(i) = ures(:,i)'*diag(whts)*ures(:,i);
        ss(i) = sqrt(1 - cs(i)^2);
        uy(:,i) = cs(i)*ux(:,i) + ss(i)*ures(:,i);
    end

%%%        checkorth(ux,ures,whts,m,k);
%%%        checkorth(ures,ures,whts,m,k);
        end
%
%
%
%
%
        function ures = whtd_res(u,whts,m,k)
%
        ures = randn(m,k);

        for j=1:k
%
        na = 2*k + 2*(j-1);

%
%        build the matrix with correct column space
%
        aa = zeros(m,na);
        ia=1;
        for i=1:j-1
%
        aa(:,ia) = ures(:,i);
        ia = ia+1;
    end

        for i=1:j-1
%
        aa(:,ia) = diag(whts) * ures(:,i);
        ia=ia+1;
    end

        for i=1:k
%
        aa(:,ia) = u(:,i);
        ia=ia+1;
    end

        for i=1:k
%
        aa(:,ia) = diag(whts)*u(:,i);
        ia=ia+1;
    end

%%%        prin2('aa=',aa(:,2),m);
%%%        prin2('ures=',ures,m);

%
%        subtract off projection of ures(:,j) from preceding
%
        kk=min(m,na);
        [ua,sa,va] = whtd_svdsmart(aa,m,na,kk);

%%%        ddd = ua*diag(sa)*va' - aa;
%%%        chk0 = norm(ddd,'fro');
%%%        prin2('chk0=',chk0,1);
        pp = ua * ua' * ures(:,j);
        ures(:,j) = ures(:,j) - pp;

        dd=whtd_scapro(ures(:,j),ures(:,j),m);
        ures(:,j) = ures(:,j) / sqrt(dd);
    end


        end
%
%
%
%
%
        function u = whtd_gramschmidt(u,m,k)
%
        u(:,1) = u(:,1) / norm(u(:,1));

        for i=2:k

        for ijk=1:2
%
        for j=1:i-1
%
        pij = sum(u(:,i) .* u(:,j));
        u(:,i) = u(:,i) - pij*u(:,j);
    end
        u(:,i) = u(:,i) / norm(u(:,i));
    end
    end

        end
%
%
%
%
%

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
