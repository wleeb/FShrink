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
