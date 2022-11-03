        function [ts,eout,einn,coutw,cinnw,cout,cinn] = ...
            whtd_emp2cheat(sy,dout,dinn,k,amus,bmus,gam,sig)
%
        ells=zeros(k,1);
        cout=zeros(k,1);
        cinn=zeros(k,1);
        ts=zeros(k,1);
%
%        parameters for each individual component
%
        for i=1:k
%
        [ts(i),cout(i),cinn(i)] = whtd_emp2pop(sy(i),gam,sig);
    end
%
%        build E matrices of weighted inner products
%
        eout = diag(1./cout) * dout * diag(1./cout);
        einn = diag(1./cinn) * dinn * diag(1./cinn);
        for i=1:k
%
        eout(i,i) = (dout(i,i) - (1-cout(i)^2)*amus(i)) / cout(i)^2;
        einn(i,i) = (dinn(i,i) - (1-cinn(i)^2)*bmus(i)) / cinn(i)^2;
    end
%
%        build C matrices of weighted inner products
%
        coutw = diag(cout) * eout;
        cinnw = diag(cinn) * einn;

        end
%
