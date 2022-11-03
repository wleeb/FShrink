        function [xhat,xfhat,errhat,khat,uyfw,tvals,etas,taus,sx,couts,...
            cinns] = fshr_core(uyfw,uyf,syf,vyf,gam0,dmu,m,n,k,iorth)
%
%        . . . memory allocation
%
        xfhat = zeros(m,n);
        xhat = zeros(m,n);
%
        taus=zeros(k,1);
        unorms=zeros(k,1);

        sx=zeros(k,1);
        cws=zeros(k,1);
        ccs=zeros(k,1);
        cinns=zeros(k,1);
        couts=zeros(k,1);
%
        uptwos = zeros(k,1);
        upones = zeros(k,1);
        dupones = zeros(k,1);
%
        etas = zeros(k,1);
        tvals = zeros(k,1);
        errs = zeros(k,1);

        errhat=0;

        gam=m/n;


%%%        prin2('dmu=',dmu,1);
%%%        prin2('gam=',gam,1);

%
%        rank estimation
%
        tol=sqrt(eps);
        khat = fshr_estrank(syf,k,gam,gam0,tol);

%%%        prinf('khat=',khat,1);
%
%        exit if rank=0
%
        if (khat == 0)
%
        return;
    end

%
%        evaluate shrinker and parameters
%
        for i=1:khat
%
        unorms(i)=norm(uyfw(:,i));
        uyfw(:,i) = uyfw(:,i) / unorms(i);
%
        [tvals(i),etas(i),errs(i),sx(i),cinns(i),couts(i),ccs(i),cws(i),upones(i),...
            dupones(i),uptwos(i),taus(i)] = fshr_emp2pars(syf(i),unorms(i),dmu,gam,gam0);

    end

%
%%%        prinstop;

        if (iorth == 1)
%
        uyfw = fshr_gramschmidt(uyfw,m,khat);
    end


%
%       reconstruct matrices
%
        xfhat = fshr_svd2mat(uyf,vyf,etas,m,n,khat);
        xhat = fshr_svd2mat(uyfw,vyf,tvals,m,n,khat);

%
%       error estimate
%
        errhat=0;
        for i=1:khat
%
        errhat = errhat + errs(i);
    end

        errhat=sqrt(errhat);

%%%        prin2('errhat=',errhat,1);



        end
