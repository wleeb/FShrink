        function ures = whtd_res(u,whts,m,k)
%
        ures = randn(m,k);

        for j=1:k
%
        na = 2*k + 2*(j-1);

%
%        build the matrix with correct column space
%
        aa = zeros(m,na);
        ia=1;
        for i=1:j-1
%
        aa(:,ia) = ures(:,i);
        ia = ia+1;
    end

        for i=1:j-1
%
        aa(:,ia) = diag(whts) * ures(:,i);
        ia=ia+1;
    end

        for i=1:k
%
        aa(:,ia) = u(:,i);
        ia=ia+1;
    end

        for i=1:k
%
        aa(:,ia) = diag(whts)*u(:,i);
        ia=ia+1;
    end

%%%        prin2('aa=',aa(:,2),m);
%%%        prin2('ures=',ures,m);

%
%        subtract off projection of ures(:,j) from preceding
%
        kk=min(m,na);
        [ua,sa,va] = whtd_svdsmart(aa,m,na,kk);

%%%        ddd = ua*diag(sa)*va' - aa;
%%%        chk0 = norm(ddd,'fro');
%%%        prin2('chk0=',chk0,1);
        pp = ua * ua' * ures(:,j);
        ures(:,j) = ures(:,j) - pp;

        dd=whtd_scapro(ures(:,j),ures(:,j),m);
        ures(:,j) = ures(:,j) / sqrt(dd);
    end


        end
%
