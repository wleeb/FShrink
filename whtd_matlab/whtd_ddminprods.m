        function tmin = whtd_ddminprods(dout,dinn,cout,cinn,tx,m,n,k,tol)
%
        aa = dout.*dinn;
        bb = cout.*cinn;

        [ainv,khat] = whtd_pseudoinv(aa,k,k,tol);
        tmin = ainv * bb * tx;

        end
%
