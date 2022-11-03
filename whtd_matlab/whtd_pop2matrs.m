        function [sy,dout,dinn,coutw,cinnw,cout,cinn] = ...
            whtd_pop2matrs(tx,eout,einn,k,amu,bmu,gam,sig);
%
        sy=zeros(1,k);
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
        dout(i,i) = cout(i)^2*eout(i,i) +  (1-cout(i)^2)*amu;
        dinn(i,i) = cinn(i)^2*einn(i,i) +  (1-cinn(i)^2)*bmu;
    end
%
%        build C matrices of weighted inner products
%
        coutw = diag(cout) * eout;
        cinnw = diag(cinn) * einn;

        end
%
