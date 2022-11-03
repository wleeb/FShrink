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
