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
