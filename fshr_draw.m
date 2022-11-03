        function [y,yf,xf,x,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,...
            ep,ep0,cov0,epw,epw0] = fshr_draw(sx,vars,m,n,k,n0,dd)
%
%        signal singular vectors
%
        ux = randn(m,k);
        ux=ux.*dd;
        ux = fshr_gramschmidt(ux,m,k);
%
        vx = randn(n,k) / sqrt(n);
%%%        vx = fshr_gramschmidt(vx,n,k);

        prin2('ux=',ux,10);
        prinf('k=',k,1)
%
%        Gaussian noise
%
        epw = randn(m,n);
        epw = epw / sqrt(n);
%%%        ep = diag(sqrt(vars)) * epw;
        ep = fshr_dleft(sqrt(vars),epw,m,n);
%
%        generate signal and spiked matrices
%
        x = ux * diag(sx) * vx';
        y = x + ep;

%
%        independent noise stream
%
        epw0 = randn(m,n0);
        epw0 = epw0 / sqrt(n0);
%%%        ep0 = diag(sqrt(vars)) * epw0;
        ep0 = fshr_dleft(sqrt(vars),epw0,m,n);

%
%        spiked F-matrix
%
        cov0 = ep0*ep0';

        cinv0=inv(cov0);
        wmat0=cinv0^(1/2);
        whts=cov0^(1/2);

        xf = wmat0 * x;
        epf = wmat0 * ep;
        yf = wmat0 * y;

        [uyf,syf,vyf] = fshr_svds(yf,m,n,k);
        [uxf,sxf,vxf] = fshr_svds(xf,m,n,k);



        if (0)
        evals = svd(ep);
        evals = evals.^2;

        prin2('evals=',evals,10);

        condn = evals(1) / evals(m);
        prin2('condn=',condn,1);
        prinstop;

    end

        evals0 = svd(ep0);
        evals0 = evals0.^2;

        prin2('evals0=',evals0,10);

        condn0 = evals0(1) / evals0(m);
        prin2('condn0=',condn0,1);

%%%        prinstop;

        end
