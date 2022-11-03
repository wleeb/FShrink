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
