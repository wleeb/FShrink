        function [ux,uy,ures,amus] = whtd_vecs(cs,m,k,whts,gam,sig)
%
        ss = zeros(k,1);
        amus = zeros(k,1);
        uy = zeros(m,k);
%
%%%        prinf('m=',m,1);
%%%        prinf('k=',k,1);

%
%        random orthonormal X vectors
%
        ux=randn(m,k);
        ux = whtd_gramschmidt(ux,m,k);
%
%        residual vectors
%
        ures = whtd_res(ux,whts,m,k);
%
%        contruct vectors for Y
%
        for i=1:k
%
        amus(i) = ures(:,i)'*diag(whts)*ures(:,i);
        ss(i) = sqrt(1 - cs(i)^2);
        uy(:,i) = cs(i)*ux(:,i) + ss(i)*ures(:,i);
    end

%%%        checkorth(ux,ures,whts,m,k);
%%%        checkorth(ures,ures,whts,m,k);
        end
%
