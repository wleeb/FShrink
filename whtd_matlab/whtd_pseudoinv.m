        function [winv,khat] = whtd_pseudoinv(w,m,n,tol)
%
        khat=min(m,n);
        [wl,ws,wr] = whtd_svdsmart(w,m,n,khat);
%%%        chk0 = norm(wl * diag(ws) * wr' - w,'fro')

        for i=1:min(m,n)
%
        rri = ws(i) / ws(1);
        if (rri < tol)
%
        khat = i - 1;
        break;
    end
    end

        wr = wr(:,1:khat);
        wl = wl(:,1:khat);
        ws = ws(1:khat);
        winv = wr * diag(1./ws) * wl';

        end
%
