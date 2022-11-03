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
