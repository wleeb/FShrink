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
