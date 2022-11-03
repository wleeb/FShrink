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
