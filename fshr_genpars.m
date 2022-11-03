        function [m,n,n0,k,rcond,sx,vars,dds,gam,gam0,taus,xcut,ycut] = ...
            fshr_genpars()
%
%        generates example problem parameters
%
        rcond = 212;

        prin2('rcond=',rcond,1);

        m=500;
        n=floor(3*m);
        n0=floor(2*m);

        k=4;
%
        prinf('m=',m,1);
        prinf('n=',n,1);
        prinf('k=',k,1);
%
        sx=zeros(k,1);
        vars = zeros(m,1);
        dds = zeros(m,k);
%
        gam=m/n;

        for i=1:m
%
        vars(i) = rcond*(m-i)/(m-1) + i/(m-1);
%%%        vars(i) = (i-1)/(m-1) + rcond*(m-i)/(m-1);
    end

%%%        vars = rcond*ones(m,1);

        vars = vars / rcond;

        prin2('vars=',vars,10);
        prin2('vars(m)=',vars(m),1);


%%%        prinstop;



%
%        signal vector weights
%
        dds=ones(m,k);

        for i=1:m
%
        dds(i,1) = sqrt(i);
        dds(i,2) = m-i+1;
    end

        for j=1:k
%
        dnorm = norm(dds(:,j));
        dds(:,j) = dds(:,j) / dnorm * sqrt(m);
    end

        prin2('dd1=',dds(:,1),10);

        chk0 = abs(norm(dds(:,1)) - sqrt(m));
        prin2('chk0=',chk0,1);


        taus = fshr_taus_exact(dds,vars,m,k);
        prin2('taus=',taus,k);


%%%        prinstop;

%%%        prinstop;

%
%        signal singular values
%
        gam0 = m/n0;
        [ycut,xcut] = fshr_cuts(gam,gam0);
        prin2('xcut=',xcut,1);


        prin2('ycut=',ycut,1);



        for i=1:k
%
        sx(i) = xcut*(1.24 + (k-i+1)/k) / sqrt(taus(i));

        sx(i) = xcut*(1.2 + (k-i+1)/k) / sqrt(taus(i));


    end

        sx(k) = 0;

%%%        sx(k)=xcut/sqrt(taus(k)) + .00001;


%%%        sx=zeros(k,1);

        prin2('sx=',sx,k);
        prin2('threshold=',xcut / sqrt(taus(1)),1);



        end
