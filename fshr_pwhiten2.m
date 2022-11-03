        function [uyfw,uyf,syf,vyf,dmu] = fshr_pwhiten2(y,ep0,n0,m,n,k)

%
%        SVD the out-of-sample noise
%
        [u0,s0,v0] = fshr_svds(ep0,m,n0,m);

%%%        chk0 = norm(ep0 - u0*diag(s0)*v0');
%%%        prin2('chk0=',chk0,1);

%
%
%        pseudo-whiten the data and perform the SVD
%
        aa = u0'*y;
        bb = fshr_dleft(1./s0,aa,m,n);
        yf = u0*bb;

        [uyf,syf,vyf] = fshr_svds(yf,m,n,k);

%
%        recolor the singular vectors (but don't normalize)
%
        aa = u0'*uyf;
        bb = fshr_dleft(s0,aa,m,n);
        uyfw = u0*bb;


%%%        prin2('yf, inside=',yf,1);

        dmu=0;
        for i=1:m
%
        dmu = dmu + s0(i)^2 / m;
    end



%%%        chk0 = norm(uyf*diag(syf) - yf*vyf,'fro') / norm(syf,'fro');
%%%        prin2('chk0=',chk0,1);

%%%        prinstop;


        end
