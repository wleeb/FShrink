        function [uyfw,uyf,syf,vyf] = fshr_pwhiten(y,cov0,m,n,k)
%
%        pseudo-whiten the data and perform the SVD
%
        whts=cov0^(1/2);

        cinv0=inv(cov0);
        wmat0=cinv0^(1/2);

%%%        chk0 = norm(wmat0*wmat0 - cinv0);
%%%        prin2('chk0=',chk0,1);

%%%        chk0 = norm(wmat0*whts - eye(m));
%%%        prin2('chk0=',chk0,1);
%%%        prinstop;

        yf = wmat0 * y;

        [uyf,syf,vyf] = fshr_svds(yf,m,n,k);
        uyfw = whts*uyf;

%%%        prin2('yf, inside=',yf,1);

%%%        chk0 = norm(uyf*diag(syf) - yf*vyf,'fro') / norm(syf,'fro');
%%%        prin2('chk0=',chk0,1);

%%%        prinstop;


        end
