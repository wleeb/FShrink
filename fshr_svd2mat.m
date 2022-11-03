        function a = fshr_svd2mat(u,v,s,m,n,k)
%
        a2 = fshr_dleft(s(1:k),v(:,1:k)',n,k);
        a = u(:,1:k) * a2;


        return;

        b = u(:,1:k) * diag(s(1:k)) * v(:,1:k)';
        chk0 = norm(a-b,'fro');
        prin2('chk0=',chk0,1);

        prinstop;

        end
