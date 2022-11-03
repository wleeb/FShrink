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
