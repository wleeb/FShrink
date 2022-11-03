        function main
        prini(13,1);


        test_basic();

%%%        compare444();


        prinstop;
        end
%
%
%
%
%
        function test_basic()
%
        randn(1,151);
        randn(1,151);
        randn(1,151);
        randn(1,151);

        [m,n,n0,k,rcond,sx,vars,dds,gam,gam0,taus,xcut,ycut] = fshr_genpars();

        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
            fshr_draw(sx,vars,m,n,k,n0,dds);


%%%        ifig=1;
%%%        fshr_wacompare(yf,gam,gam0,ifig);


%%%        prinstop;

        prin2('vars=',vars,10);
        prin2('sx=',sx,k);
%
%        denoise with pseudo-whitening
%
        iorth=0;

        [xhat,xfhat,errhat,khat,uyfw,vy,tvals,etas,taus_hat,sx_hat,couts,...
            cinns,dmu] = fshr_approx2(y,ep0,n0,m,n,k,iorth);


        prin2('xhat=',xhat,100);
        prin2('xfhat=',xfhat,100);
        prin2('tval=',tvals,khat);
        prin2('taus_hat=',taus_hat,khat);
        prin2('cinns=',cinns,khat);
        prin2('couts=',couts,khat);
        prin2('uyfw=',uyfw,100);
        prin2('vy=',vy,100);


        prinf('khat=',khat,1);


        err = norm(xhat - xmat,'fro');

        prin2('errhat=',errhat,1);
        prin2('error=',err,1);


        if (1)

        [xhat,xfhat,errhat,khat,uyfw,vy,tvals,etas,taus_hat,sx_hat,couts,...
            cinns,dmu] = fshr_approx(y,cov0,gam0,m,n,k,iorth);


        prin2('xhat=',xhat,100);
        prin2('xfhat=',xfhat,100);
        prin2('tval=',tvals,khat);
        prin2('taus_hat=',taus_hat,khat);
        prin2('cinns=',cinns,khat);
        prin2('couts=',couts,khat);
        prin2('uyfw=',uyfw,100);
        prin2('vy=',vy,100);


        prinf('khat=',khat,1);


        err = norm(xhat - xmat,'fro');

        prin2('errhat=',errhat,1);
        prin2('error=',err,1);

    end

%
%        compare to exact least squares solution
%
        [xhat2,xfhat2,errhat2,tvals2,etas2] = fshr_exact(y,xmat,...
            cov0,whts,m,n,n0,k);


        prin2('taus_hat=',taus_hat,k);
        prin2('taus, true=',taus,k);

        prin2('errhat=',errhat,1);
        prin2('errhat2=',errhat2,1);

        prin2('etas=',etas,k);
        prin2('topt2=',etas2,k);

        prin2('tvals=',tvals,k);
        prin2('tvals2=',tvals2,k);


        prinstop;
        end
%
%
%
%
%
        function compare444()

        fname = 'compare_cosines';


        gams0 = [.1,.3,.5];
        ngams0 = length(gams0);

        gam=.5;

        m=2000;


        rcond_max=100;

%
%        set parameters, for fixed beta
%
        nconds = 100;

        rconds = zeros(nconds,1);
        couts = zeros(nconds,ngams0);
        couts_or = zeros(nconds,1);
        couts_nt = zeros(nconds,1);


        for i=1:nconds
%
        rconds(i) = 1 + (rcond_max - 1)*(i-1)/(nconds-1);

    end


        prin2('rconds=',rconds,nconds);


%%%        prinstop;
%
%        detection threshold for all methods
%
        gmax0 = max(gams0);

        prin2('gmax0=',gmax0,1);

        [ycut,xcut] = fshr_cuts(gam,gmax0);

        prin2('xcut=',xcut,1);
        prin2('ycut=',ycut,1);

        vars_max = make_vars(rcond_max,m);
        xcut_bbp = eval_bbp(vars_max,m,gam);

        prin2('bbp=',xcut_bbp,1);
        prin2('xcut=',xcut,1);


%
%        signal singular value
%
        sx = max(xcut_bbp,xcut) + .5;

        prin2('sx=',sx,1);


        prin2('vars=',vars_max,10);
        prin2('vars(m)=',vars_max(m),1);


%
%        asymptotic inner products
%

        for i=1:nconds
%

        vars = make_vars(rconds(i),m);

        [couts(i,1),couts_or(i),couts_nt(i)] = cosines_all(sx,vars,m,gam,gams0(1));



        for j=1:ngams0
%
        dds = ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);
        dmu = sum(vars) / m;
%
        [tval,eta,err,syf,cinn,couts(i,j),cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gams0(j),tau);
    end

    end



        prinar2('couts=',couts,nconds,ngams0);
        prin2('couts, oracle=',couts_or,nconds);
        prin2('couts, no=',couts_nt,nconds);



        save(fname,'couts','couts_or','couts_nt','m','gam','ngams0',...
            'gams0','nconds','rconds','sx','xcut','xcut_bbp');


        plot_cosines555(' ',fname);

        end
%
%
%
%
%
        function plot_cosines555(ptitle,fname)

        load(fname);
        who;

        fh = figure();
        hold on;
        box on;
        grid on;

        legs = strings(ngams0+2,1);

        legs(1) = '$\beta=0$ (oracle whitening)'


        plot(rconds,couts_or,'b-.','linewidth',2);

        for j=1:ngams0
%
        plot(rconds,couts(:,j),'linewidth',2);

        sss = int2str(10*gams0(j));

        legs(j+1) = strcat('$\beta=0.',sss,'$');
    end


        legs(ngams0+2) = 'No whitening';


        plot(rconds,couts_nt,'r--','linewidth',2);

%%%        plot(log2(rconds),cmeans-cmeans_or,'ks-','linewidth',2);
%%%        plot(log2(rconds),cmeans_nt-cmeans_or,'ro--','linewidth',2);


        xlim([rconds(1),rconds(nconds)]);


        legend(legs,'interpreter','latex','fontsize',20,'location','southwest');
        set(fh,'Position',[500,500,900,800]);


        xlabel('Condition number','interpreter','latex','fontsize',20);
        ylabel('Cosine','interpreter','latex','fontsize',20);

        title(ptitle,'interpreter','latex','fontsize',20)

        set(fh,'PaperPositionMode','auto')
        print(fh,fname,'-depsc','-r0')





        end
%
%
%
%
%
        function [cout,cout_or,cout_nt] = cosines_all(sx,vars,m,gam,gam0)


%
%        pseudo-whitened angle
%
        dds = ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);
        dmu = sum(vars) / m;
%
        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gam0,tau);


%
%        oracle whitened
%
        [tval,eta,err,syf,cinn,cout_or,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,0,tau);


%
%        no whitening
%
        ell=sx^2;
        awhts=ones(m,1)/m;
        n=1;
        bs=1;
        bwhts=1;
        [rlam,cout_nt,cinn] = mpbdry_sforw(ell,vars,bs,awhts,bwhts,m,n,gam);


        end
%
%
%
%
%
        function [cprods,cprods_or,cprods_nt] = angle_sweep(rconds,nconds,sx,m,n,n0)

        prin2('rconds=',rconds,nconds);
%%%        prinstop;


        cprods = zeros(nconds,1);
        cprods_or = zeros(nconds,1);
        cprods_nt = zeros(nconds,1);

%
%        set random parameters
%
        [xmat,ux,vx,epw,epw0] = draw_enough(sx,m,n,n0);


        for ijk=1:nconds
%

        vars = make_vars(rconds(ijk),m);

%%%        vars = vars / vars(1);

        prin2('vars=',vars,10);
        prin2('vars(m)=',vars(m),1);

%
%        find all inner products
%
        [cprod,cprod_or,cprod_nt] = angles_all(xmat,epw,epw0,vars,ux,vx,m,n,n0);


        cprods(ijk) = cprod;
        cprods_or(ijk) = cprod_or;
        cprods_nt(ijk) = cprod_nt;


        prin2('cprod=',cprod,1);
        prin2('cprod_or=',cprod_or,1);
        prin2('cprod_nt=',cprod_nt,1);

    end

        prin2('cprods=',cprods,nconds);
        prin2('cprods_or=',cprods_or,nconds);
        prin2('cprods_nt=',cprods_nt,nconds);



        end
%
%
%
%
%
        function plot_cosines(ptitle,fname)

        load(fname);
        who;

        fh = figure();
        hold on;
        box on;
        grid on;

        plot(log2(rconds),cmeans_or,'b<-.','linewidth',2);
        plot(log2(rconds),cmeans,'ks-','linewidth',2);
        plot(log2(rconds),cmeans_nt,'ro--','linewidth',2);


%%%        plot(log2(rconds),cmeans-cmeans_or,'ks-','linewidth',2);
%%%        plot(log2(rconds),cmeans_nt-cmeans_or,'ro--','linewidth',2);


        xlim([log2(rconds(1)),log2(rconds(nconds))]);


        legend('Oracle whitening','Pseudo-whitening','No whitening',...
            'interpreter','latex','fontsize',20,'location','southwest');
        set(fh,'Position',[500,500,900,800]);


        xlabel('$\log_2$ condition number','interpreter','latex','fontsize',20);
        ylabel('Cosine','interpreter','latex','fontsize',20);

        title(ptitle,'interpreter','latex','fontsize',20)

        set(fh,'PaperPositionMode','auto')
        print(fh,fname,'-depsc','-r0')





        end
%
%
%
%
%
        function compare_angles()
%
        iplot=0;


        gams0=zeros(10,1);
        ptitles=strings(10,1);
        fnames=strings(10,1);


        gams0=[1/2, 3/4, 1/4, 1/8];

        i=1;
        ptitles(i) = '$\beta$ = 1/2';
        fnames(i) = 'angles_half';

        i=i+1;
        ptitles(i) = '$\beta$ = 3/4';
        fnames(i) = 'angles_threefourths';

        i=i+1;
        ptitles(i) = '$\beta$ = 1/4';
        fnames(i) = 'angles_quarter';

        i=i+1;
        ptitles(i) = '$\beta$ = 1/8';
        fnames(i) = 'angles_eighth';



        if (iplot)
%

        for i=1:4


        plot_cosines(ptitles(i),fnames(i));


    end

        prinstop;
    end


%

%

        ptitles
        fnames


        for i=1:4
%
        compare_angles0(gams0(i),fnames(i));

    end
%%%        plot_cosines(ptitle,fname);
        prinstop;





        prinstop;
%%%        test_maxbeta444();
        plot_bbp_rconds();

%%%        plot_bbp_alphas();

        prinstop;
        end
%
%
%
%
%
        function compare_angles0(gam0,fname)

        fname = strcat(fname,'.mat')
%
%        set parameters, for fixed beta
%
        nconds = 12;
        ndraws = 3;

        rconds = zeros(nconds,1);
        cprods = zeros(nconds,ndraws);
        cprods_or = zeros(nconds,ndraws);
        cprods_nt = zeros(nconds,ndraws);


        for i=1:nconds
%
%%%        rconds(i) = 1 + (rcond_max - 1)*(i-1)/(nconds-1);
        rconds(i) = 2^(i-1);

    end



        prin2('rconds=',rconds,nconds);

%%%        prinstop;


        gam=1/4;


        m=600;
        n=floor(m/gam);
        n0=floor(m/gam0);

        prinf('m=',m,1);
        prinf('n=',n,1);
        prinf('n0=',n0,1);
%

%%%        prinstop;
%
%        signal singular values
%
        [ycut,xcut] = fshr_cuts(gam,gam0);

        prin2('xcut=',xcut,1);
        prin2('ycut=',ycut,1);

%%%        sx = xcut*1.5;

        sx = xcut + 2;

        prin2('sx=',sx,1);
        prin2('threshold=',xcut,1);


        for iii=1:ndraws

        prinf('iii=',iii,1);

        rng(iii);

        [cprods(:,iii),cprods_or(:,iii),cprods_nt(:,iii)] = ...
            angle_sweep(rconds,nconds,sx,m,n,n0);
    end

        prinar2('cprods=',cprods,nconds,3);
        prinar2('cprods=',cprods_or,nconds,3);
        prinar2('cprods=',cprods_nt,nconds,3);


        cmeans = mean(cprods,2);
        cmeans_or = mean(cprods_or,2);
        cmeans_nt = mean(cprods_nt,2);

        prin2('cmeans=',cmeans,nconds);
        prin2('cmeans_or=',cmeans_or,nconds);
        prin2('cmeans_nt=',cmeans_nt,nconds);

%%%        prinar2('cprods=',cprods,nconds,ndraws);
%%%        prinar2('cprods=',cprods_or,nconds,ndraws);
%%%        prinar2('cprods=',cprods_nt,nconds,ndraws);


        save(fname,'cmeans','cmeans_or','cmeans_nt','m','n','n0','gam',...
            'gam0','ndraws','nconds','rconds','sx','xcut');



%%%        plot_cosines();

        end
%
%
%
%
%
        function [x,ux,vx,epw,epw0] = draw_enough(sx,m,n,n0)
%
        k=1;
        dd=ones(m,1);

%
%        signal singular vectors
%
        ux = randn(m,1);
        ux = ux / norm(ux);
%
        vx = randn(n,k) / sqrt(n);

        prin2('ux=',ux,10);
        prinf('k=',k,1)
%
%        Gaussian noise
%
        epw = randn(m,n);
        epw = epw / sqrt(n);
%
        x = ux * diag(sx) * vx';
%
%        independent noise stream
%
        epw0 = randn(m,n0);
        epw0 = epw0 / sqrt(n0);

        end
%
%
%
%
%
        function [cprod,cprod_or,cprod_nt] = angles_all(x,epw,epw0,vars,ux,vx,m,n,n0)
%
%        returns PC angles for all three methods (pseudo-whitening + recoloring,
%        oracle whitening + recoloring, and no transformation)
%
        k=1;

%
%        color in-sample noise and build matrix
%
        ep = fshr_dleft(sqrt(vars),epw,m,n);
        y = x + ep;

%
%        color out-of-sample noise
%
        ep0 = fshr_dleft(sqrt(vars),epw0,m,n);
        cov0 = ep0*ep0';

        cinv0=inv(cov0);
        wmat0=cinv0^(1/2);
        whts=cov0^(1/2);

        xf = wmat0 * x;
        epf = wmat0 * ep;
        yf = wmat0 * y;

        [uyf,syf,vyf] = fshr_svds(yf,m,n,k);
        [uxf,sxf,vxf] = fshr_svds(xf,m,n,k);



%
%        recolor the singular vector
%
        uhat = whts*uyf;
        uhat = uhat / norm(uhat);

        chk0 = norm(uhat) - 1;
        prin2('chk0=',chk0,1);

        cprod = abs(sum(uhat.*ux));


%
%        evaluate estimated vectors using oracle whitening
%
        [uhat_or,vhat_or] = color_oracle(y,x,ep,epw,vars,m,n);
        cprod_or = abs(sum(uhat_or.*ux));



%
%        inner products with no transformation
%
        [uy,sy,vy] = fshr_svds(y,m,n,k);
        cprod_nt = abs(sum(uy.*ux));

        chk0 = norm(y*vy - uy*diag(sy));


        prin2('cprod=',cprod,1);
        prin2('cprod_or=',cprod_or,1);






        if (0)

%
%        check projection
%
        gam0=m/n0;

        iorth=0;
        [xhat,xfhat,errhat,khat,uyfw,vyf,tvals,etas,taus,sx,couts,...
            cinns,dmu] = fshr_approx(y,cov0,gam0,m,n,k,iorth);

        chk0 = norm(xhat - uhat*uhat'*xhat);
        prin2('chk0=',chk0,1);

%
%        evaluate cosines
%

        prin2('cprod=',cprod,1);
        prin2('couts=',couts,1);
        dif = cprod - couts;
        prin2('dif=',dif,1);

    end





        if (0)
%
%        check projection
%
        [xhat4,khat4] = oracle_wshr(y,m,n,k,vars);
        chk0 = norm(xhat4 - uhat_or*uhat_or'*xhat4);
        prin2('chk0=',chk0,1);


        [uy2,sy2,vy2] = fshr_svds(y,m,n,m);
        chk0 = norm(y*vy2 - uy2*diag(sy2));

        chk0 = 1 - abs(sum(uy.*uy2(:,1)));
        chk0 = chk0 + 1 - abs(sum(vy.*vy2(:,1)));
        prin2('chk0=',chk0,1);
    end




        end
%
%
%
%
%
        function [uhat,vhat] = color_oracle(y,x,ep,epw,vars,m,n)
%
%        apply oracle whitening and recoloring to top vector
%
        yw = diag(1./sqrt(vars))*y;


        if (0)
        xw = diag(1./sqrt(vars))*x;

        chk0 = norm(yw - xw - epw);
        prin2('chk0=',chk0,1);

        chk0 = norm(y - x - ep);
        prin2('chk0=',chk0,1);

    end

%
%        decompose whitened matrix
%
        k=1;
        [uyw,syw,vhat] = fshr_svds(yw,m,n,k);

        chk0 = norm(yw*vhat - uyw*syw,'fro');
        prin2('chk0=',chk0,1);

%
%        rewhiten singular vector
%
        uhat = diag(sqrt(vars)) * uyw;
        unorm = norm(uhat);

        uhat = uhat / unorm;

        end
%
%
%
%
%
        function plot_bbp_alphas()

        nalps = 200;
        alp_max = 8;
        alp_min = 0;

        ngams = 5;

        bbps = zeros(nalps,ngams);
        gopts0 = zeros(nalps,ngams);


        gams = zeros(ngams,1);

        gams(1)=100;
        gams(2)=10;
        gams(3)=1;
        gams(4)=.1;
        gams(5)=.01;


        prin2('gams=',gams,ngams);
%%%        prinstop;




        gam=.01;
        m=2000;
        vars=make_vars2(5*alp_max/nalps,m);
        [gopt0,xcut_bbp,tau] = betathresh_secant(vars,m,gam);

        prin2('gopt0=',gopt0,1);

        [gopt0_2,xcut_bbp_2,tau_2] = betathresh_fmla(vars,m,gam);

        prin2('gopt0',gopt0,1);
        prin2('gopt0_2',gopt0_2,1);

        chk0 = gopt0 - gopt0_2;
        prin2('chk0=',chk0,1);

%%%        prinstop;


        for i=1:ngams
%
        prinf('i=',i,1);
        prin2('gamma=',gams(i),1);

%%%        prini(0,0);
        [gopts0(:,i),bbps(:,i),alps] = plot_bbp_alphas0(gams(i),nalps,alp_max,alp_min);
        prini(13,0);
    end



        prin2('gopts0=',gopts0,10);
        prin2('bbps=',bbps,10);


%
%        plot bbps
%
        figure(1);
        hold on;
        box on;
        grid on;

        for i=1:ngams
        plot(alps,bbps(:,i),'linewidth',2);

        sss = int2str(log10(gams(i)));
        legs(i) = strcat("$\gamma=10^{",sss,"}$");
    end

        legend(legs,'interpreter','latex','fontsize',20,'location','northwest');
%%%        legend(legs,'interpreter','latex','fontsize',20,'location','eastoutside');
        set(figure(1),'Position',[500,500,900,800]);

        xlim([alp_min,alp_max]);
        ylim([0,max(bbps(:,1))+1.3]);


        xlabel('$\alpha$','interpreter','latex','fontsize',20);
        ylabel('BBP','interpreter','latex','fontsize',20);



        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'bbp_alphas','-depsc','-r0')


%%%        prinstop;




        figure(2);
        hold on;
        box on;
        grid on;

        for i=1:ngams
%
        prinf('i=',i,1);
        plot(alps,gopts0(:,i),'linewidth',2);

        sss = int2str(log10(gams(i)));
        legs(i) = strcat("$\gamma=10^{",sss,"}$");
    end


         gopts0


        legend(legs,'interpreter','latex','fontsize',20,'location','southeast');
%%%        legend(legs,'interpreter','latex','fontsize',20,'location','eastoutside');
        set(figure(2),'Position',[500,500,900,800]);

        xlim([alp_min,alp_max]);
        ylim([0,max(gopts0(:,1))]);


        xlabel('$\alpha$','interpreter','latex','fontsize',20);
        ylabel('$\beta$','interpreter','latex','fontsize',20);



        set(figure(2),'PaperPositionMode','auto')
        print(figure(2),'betas_thresh_alphas','-depsc','-r0')






        end
%
%
%
%
%
        function [gopt0,xcut_bbp,tau] = betathresh_fmla(vars,m,gam)
%
%        bbp transition and tau
%
        xcut_bbp = eval_bbp(vars,m,gam);
%
        dds=ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);

%%%        prin2('tau=',tau,1);
%%%        prin2('xcut_bbp=',xcut_bbp,1);

        aa = (xcut_bbp^2*tau+1)^2;
        bb = -1+gam - 2*xcut_bbp^2*tau*xcut_bbp^2*tau -2*xcut_bbp^2*tau;
        cc = xcut_bbp^4*tau^2-gam;


        gopt0 = (-bb-sqrt(bb^2 - 4*aa*cc)) / (2*aa);

        end
%
%
%
%
%
        function [gopts0,bbps,alps] = plot_bbp_alphas0(gam,nalps,alp_max,alp_min)
%
%        plots BBP and maximum beta as functions of
%        alpha, for variances scaling like t^alpha
%
        alps = zeros(nalps,1);

        bbps = zeros(nalps,1);
        gopts0 = zeros(nalps,1);


        m=2000;

%
%        plot BBP as function of condition number
%

        for i=1:nalps
%
        prinf('i=',i,1);
        alps(i) = alp_min + (alp_max-alp_min)*(i-1)/(nalps-1);
    end


        prin2('alphas=',alps,nalps);


%
%        bbp and maximum beta, as function of condition number
%
        for i=1:nalps
%
        vars = make_vars2(alps(i),m);

        prin2('vars=',vars,10);

%%%        plot(vars)

%%%        prinstop;


%%%        rr=vars(m)/vars(1);
%%%        prin2('rr=',rr,1);
%%%        prin2('rconds(i)=',rconds(i),1);

        bbps(i) = eval_bbp(vars,m,gam);

        prini(0,0);
        [gopts0(i),xcut_bbp,tau] = betathresh_fmla(vars,m,gam);
        prini(13,0);

        chk0 = bbps(i) - xcut_bbp;
        prin2('chk0=',chk0,1);

%%%        prinstop;

    end


%%%        plot(vars);
%%%        prinstop;

        prin2('bbps=',bbps,nalps);

        prin2('gopts0=',gopts0,nalps);


        end
%
%
%
%
%
        function [gopt0,gams0,vals,ierr] = betathresh_dumb(vars,m,gam,nn)
%
%        . . . finds largest gamma0 (beta) where sigma(gamm0) = BBP
%

        gams0 = zeros(nn,1);
        xcuts = zeros(nn,1);
        vals = zeros(nn,1);

        ierr=0;
%
%        bbp transition
%
        xcut_bbp = eval_bbp(vars,m,gam);

        gmin0 = 0;
        gmax0 = .99;

        prin2('gmin0=',gmin0,1);
        prin2('gmax0=',gmax0,1);

        dds=ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);

%
%        evaluate cutoff as function of gamma0
%
        for i=1:nn
%

        prinf('i=',i,1);

        gams0(i) = gmin0 + (gmax0-gmin0)*(i-1)/(nn-1);

        xcuts(i) = fshr_xcut_tau(gam,gams0(i),tau);
        vals(i) = xcuts(i) - xcut_bbp;

    end

        if (vals(nn) <= 0)
%
        gopt0=1;
        return;
    end


%%%        plot(gams0,vals);


        gopt0=-1;

%
%        brute force search for root
%
        for i=1:nn-1
%
        if (vals(i) <= 0 && vals(i+1) >= 0)
%
        gopt0 = (gams0(i) + gams0(i+1)) / 2;
        break;
        end
    end

        prin2('gopt0=',gopt0,1);

        if (gopt0 < 0)
%
        ierr=1;
    end


        end
%
%
%
%
%
        function plot_bbp_rconds()

        ngams = 5;

        nconds = 200;
        rmax = 500;
        rmin = 1;
        gam = .5;

        bbps = zeros(nconds,ngams);
        gopts0 = zeros(nconds,ngams);


        gams = zeros(ngams,1);

        gams(1)=100;
        gams(2)=10;
        gams(3)=1;
        gams(4)=.1;
        gams(5)=.01;





        for i=1:ngams
%
        prinf('i=',i,1);
        prin2('gamma=',gams(i),1);

        prini(0,0);
        [bbps(:,i),gopts0(:,i),rconds] = plot_bbp_rconds0(nconds,rmax,rmin,gams(i));
        prini(13,0);
    end

%%%        bbps
%%%        gopts0

%%%        prinstop;

        prin2('gopts0=',gopts0,10);
        prin2('bbps=',bbps,10);


%
%        plot bbps
%
        figure(1);
        hold on;
        box on;
        grid on;

        for i=1:ngams
        plot(rconds,bbps(:,i),'linewidth',2);

        sss = int2str(log10(gams(i)));
        legs(i) = strcat("$\gamma=10^{",sss,"}$");
    end

        legend(legs,'interpreter','latex','fontsize',20,'location','northwest');
%%%        legend(legs,'interpreter','latex','fontsize',20,'location','eastoutside');
        set(figure(1),'Position',[500,500,900,800]);

        xlim([rmin,rmax]);
        ylim([0,max(bbps(:,1))+1.3]);


        xlabel('$\kappa$','interpreter','latex','fontsize',20);
        ylabel('BBP','interpreter','latex','fontsize',20);



        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'bbp_conds','-depsc','-r0')



%%%        prinstop;




        figure(2);
        hold on;
        box on;
        grid on;

        for i=1:ngams
        plot(rconds,gopts0(:,i),'linewidth',2);
    end



        legend(legs,'interpreter','latex','fontsize',20,'location','southeast');
%%%        legend(legs,'interpreter','latex','fontsize',20,'location','eastoutside');
        set(figure(2),'Position',[500,500,900,800]);

        xlim([rmin,rmax]);
        ylim([0,max(gopts0(:,1))]);


        xlabel('$\kappa$','interpreter','latex','fontsize',20);
        ylabel('$\beta$','interpreter','latex','fontsize',20);



        set(figure(2),'PaperPositionMode','auto')
        print(figure(2),'betas_thresh_conds','-depsc','-r0')




        end
%
%
%
%
%
        function [bbps,gopts0,rconds] = plot_bbp_rconds0(nconds,rmax,rmin,gam)
%
%        plots BBP and maximum beta as functions of
%        the condition number, for equispaced variances
%
        rconds = zeros(nconds,1);
        bbps = zeros(nconds,1);
        gopts0 = zeros(nconds,1);


        m=2000;

%
%        plot BBP as function of condition number
%

        for i=1:nconds
%
        prinf('i=',i,1);
        rconds(i) = rmin + (rmax-rmin)*(i-1)/(nconds-1);
    end


        prin2('rconds=',rconds,nconds);

%
%        bbp and maximum beta, as function of condition number
%
        for i=1:nconds
%
        vars = make_vars(rconds(i),m);

        prin2('vars=',vars,10);

        rr=vars(m)/vars(1);
        prin2('rr=',rr,1);
        prin2('rconds(i)=',rconds(i),1);

        bbps(i) = eval_bbp(vars,m,gam);
        [gopts0(i),xcut_bbp,tau] = betathresh_fmla(vars,m,gam);

        chk0 = bbps(i) - xcut_bbp;
        prin2('chk0=',chk0,1);

    end
        prin2('bbps=',bbps,nconds);

        prin2('gopts0=',gopts0,nconds);


        end
%
%
%
%
%
        function test_maxbeta444()

        m=1000;

        k=1;

%
%        detectable signal singular value
%

        gam=.1;
        gam0=.75;

        rcond=100;

        alp=5;
        vars = make_vars2(alp,m);
        [xcut,xcut2,xcut3] = cutoff_both(vars,m,gam,gam0);

        prin2('vars=',vars,m);

        rr = vars(m) / vars(1);

        prin2('rr=',rr,1);

        prin2('xcut=',xcut,1);
        prin2('xcut2=',xcut2,1);
        prin2('xcut3=',xcut3,1);



        nn=1000;

        [gopt0_old,gams0,vals,ierr] = betathresh_dumb(vars,m,gam,nn);

        prin2('gopt0=',gopt0_old,1);

%%%        prinstop;

        [gopt0,xcut_bbp,tau] = betathresh_secant(vars,m,gam);

        prin2('gopt0=',gopt0_old,1);
        prin2('gopt0, returned=',gopt0,1);


        [val,xcut] = thresh_eval555(gam,gopt0,tau,xcut_bbp);


        prin2('val=',val,1);

%%%        return;
%%%        prinstop;

%%%        prin2('gams0=',gams0,nn);
%%%        prin2('vals=',vals,nn);



        dzero = 0;

        hold on;
        plot(gams0,vals);
        plot(gopt0,dzero,'*');

        plot(gams0,dzero*ones(nn,1),'k');



        prinstop;
        end
%
%
%
%
%
        function [gopt0,xcut_bbp,tau] = betathresh_secant(vars,m,gam)


%
%        bbp transition and tau
%
        xcut_bbp = eval_bbp(vars,m,gam);
%
        dds=ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);


%
%        if variances are constant, optimal beta is 0
%
        vmin = min(vars);
        vmax = max(vars);
        if (vmin == vmax)
%
        gopt0=0;
        prin2('inside secant code, gopt0=',gopt0,1);
        return;
    end

        niter=20;

        vals = zeros(niter,1);
        xcuts = zeros(niter,1);
        gams0 = zeros(niter,1);


        prin2('gams0=',gams0,5);


%
%        initialize the two values of beta
%
        gams0(1) = .1;

        for ijk=1:1000
%
        [vals(1),xcuts(1)] = thresh_eval555(gam,gams0(1),tau,xcut_bbp);

        if (vals(1) < 0)
%
        gams0(1) = (1+gams0(1))/2;
        prinf('ijk=',ijk,1);

        prin2('gams0(1)=',gams0(1),1);
        continue;
    end
        break;
    end

        gams0(2) = (1+gams0(1))/2;;


        prin2('gams0(1)=',gams0(1),1);
        prin2('vals(1)=',vals(1),1);


        [vals(1),xcuts(1)] = thresh_eval555(gam,gams0(1),tau,xcut_bbp);
        [vals(2),xcuts(2)] = thresh_eval555(gam,gams0(2),tau,xcut_bbp);

        prin2('vals=',vals,5);
        prin2('xcuts=',xcuts,5);
        prin2('gams0=',gams0,5);


%%%        prinstop;


        thresh = sqrt(eps);
        prin2('thresh=',thresh,1);
%%%        prinstop;

        kplus=0;
%
%        secant updates
%
        for i=3:niter

        grad = (vals(i-1) - vals(i-2)) / (gams0(i-1) - gams0(i-2));

        gams0(i) = gams0(i-1) - vals(i-1)/grad;

%%%        xcuts(i) = fshr_xcut_tau(gam,gams0(i),tau);
%%%        vals(i) = xcuts(i) - xcut_bbp;

        [vals(i),xcuts(i)] = thresh_eval555(gam,gams0(i),tau,xcut_bbp);

%
%        exit if outside interval
%
        if (gams0(i) < 0 || gams0(i) > 1)
%
        break;
    end

%
%        check convergence
%
        if (abs(vals(i)) < thresh)
%
        kplus = kplus + 1;
        prinf('incrementing, kplus=',kplus,1);
        prin2('incrementing, vals(i)=',vals(i),1);
    end

        if (kplus == 2)
%
        prinf('exiting, kplus=',kplus,1);
        break;
    end

        prin2('grad=',grad,1);
    end

        prin2('gams0=',gams0,niter);
        prin2('vals=',vals,niter);

        gopt0 = gams0(i);

        prin2('gopt0=',gopt0,1);

%%%        nn=100;
%%%        [gopt0_2,gams0,vals,ierr] = betathresh_dumb(vars,m,gam,nn);
%%%        plot(gams0,vals);

        prin2('gam=',gam,1);

        prin2('vars',vars,10);

        prinf('kplus=',kplus,1);


        [val22,xcut22] = thresh_eval555(gam,gopt0,tau,xcut_bbp);


%%%        val22
%%%        prinstop;
%%%        vals
%%%        gams0


        if (kplus == 2)
%
        prinf('converged, kplus=',kplus,1);
        return;
    end

%
%        if did not converge, approach from other side
%
        gams0(1) = 0;
        [vals(1),xcuts(1)] = thresh_eval555(gam,gams0(1),tau,xcut_bbp);


%
%        initialize to the left of root
%
        gg0 = .1;

        for i=1:100
%
        [val,xcut] = thresh_eval555(gam,gg0,tau,xcut_bbp);

        prin2('gg0=',gg0,1);

        if (val > 0)
%
        gg0 = gg0/2;
        continue
    end

        break;
    end

%%%        prinstop;
        gams0(2) = gg0;

        prin2('gams0=',gams0,3);

        [vals(2),xcuts(2)] = thresh_eval555(gam,gams0(2),tau,xcut_bbp);

        prin2('vals=',vals,5);
        prin2('xcuts=',xcuts,5);
        prin2('gams0=',gams0,5);

%%%        prinstop;


        kplus=0;
%
%        secant updates
%
        for i=3:niter

        grad = (vals(i-1) - vals(i-2)) / (gams0(i-1) - gams0(i-2));

        gams0(i) = gams0(i-1) - vals(i-1)/grad;

%%%        xcuts(i) = fshr_xcut_tau(gam,gams0(i),tau);
%%%        vals(i) = xcuts(i) - xcut_bbp;

        [vals(i),xcuts(i)] = thresh_eval555(gam,gams0(i),tau,xcut_bbp);

        if (abs(vals(i)) < thresh)
%
        kplus = kplus + 1;
    end

        if (kplus == 2)
%
        break;
    end

        prin2('grad=',grad,1);
    end


        gopt0 = gams0(i);
        [val,xcut] = thresh_eval555(gam,gopt0,tau,xcut_bbp);

        prin2('vals=',val,1);
        prin2('xcut=',xcut,1);
        prin2('gopt0=',gopt0,1);


%%%%        prinstop;



        end
%
%
%
%
%
        function xcut_bbp = eval_bbp(vars,m,gam)

        n=1;
        dzero = 0;
        awhts=ones(1,m) / m;
        bwhts=1;
        bs=1;

        prinf('n=',n,1);
        prinf('m=',m,1);
        prin2('gam=',gam,1);
%%%        prinstop;
        [ell,bedge,err] = mpbdry_thresh(vars,bs,awhts,bwhts,m,n,gam);

        xcut_bbp = sqrt(ell);


        end
%
%
%
%
%
        function [val,xcut] = thresh_eval555(gam,gam0,tau,xcut_bbp)
%

        xcut = fshr_xcut_tau(gam,gam0,tau);
        val = xcut - xcut_bbp;

        end
%
%
%
%
%
        function plot_errbetas_alphas()
%
%        evaluates maximum beta where new method beats OptShrink,
%        for polynomial variances with specified exponent
%        and varying values of gamma
%

        ngams = 5;
        nalps = 200;

        gams = zeros(ngams,1);
        gopts0 = zeros(nalps,ngams);

        gams(1) = 100;
        gams(2) = 10;
        gams(3) = 1;
        gams(4) = .1;
        gams(5) = .01;


        prin2('gams=',gams,ngams);




        alp_max = 8;
        alp_min = 0;


        for i=1:ngams

        prinf('i=',i,1);
        prin2('gamma=',gams(i),1);

%%%        prini(0,0);
        [gopts0(:,i),alps] = plot_errbetas_alphas0(gams(i),alp_max,alp_min,nalps);
        prini(13,0);
    end


        hold on;
        box on;
        grid on;

        for i=1:ngams
%
        plot(alps,gopts0(:,i),'linewidth',2);

        sss = int2str(log10(gams(i)));
        legs(i) = strcat("$\gamma=10^{",sss,"}$");
    end



        legend(legs,'interpreter','latex','fontsize',20,'location','southeast');
        set(figure(1),'Position',[500,500,900,800]);

        xlim([alp_min,alp_max]);
        ylim([min(gopts0(:,ngams)),max(gopts0(:,1))+.01]);


        xlabel('$\alpha$','interpreter','latex','fontsize',20);
        ylabel('$\beta$','interpreter','latex','fontsize',20);



        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'betas_errs_alphas','-depsc','-r0')





        prinstop;
        end
%
%
%
%
%
        function [gopts0,alps] = plot_errbetas_alphas0(gam,alp_max,alp_min,nalps)
%
%        evaluates maximum beta where new method beats OptShrink,
%        for polynomial variances with specified exponent
%


        alps = zeros(nalps,1);
        bbps = zeros(nalps,1);
        gopts0 = zeros(nalps,1);


        m=2000;


        prinf('m=',m,1);
        prin2('gam=',gam,1);


        for i=1:nalps
%
        prinf('i=',i,1);
        alps(i) = alp_min + (alp_max-alp_min)*(i-1)/(nalps-1);
    end


        prin2('alphas=',alps,nalps);

        alpmax = alps(nalps);
        vars_max = make_vars2(alpmax,m);

        prinf('m=',m,1);


%%%        prini(0,0);
%
%        evaluate error as function of condition number
%


        for i=1:nalps
%
        vars = make_vars2(alps(i),m);
        xcut_bbp = eval_bbp(vars,m,gam);

        prin2('bbp=',xcut_bbp,1);


        sx = 1.1*xcut_bbp;


        prin2('vars=',vars,10);

        [gopts0(i),err_opt,dmu,tau,ierr] = maxbeta_err_secant(sx,vars,m,gam);
    end

%%%        prini(13,0);

        prin2('gopts0=',gopts0,10);

        end
%
%
%
%
%
        function plot_errbetas_rconds()
%
%        evaluates maximum beta where new method beats OptShrink,
%        for equispaced variances with specified condition numbers
%        and varying values of gamma
%

        ngams = 5;
        nconds = 200;

%%%        nconds = 10;


        gams = zeros(ngams,1);


        gams(1) = 100;
        gams(2) = 10;
        gams(3) = 1;
        gams(4) = .1;
        gams(5) = .01;

        prin2('gams=',gams,ngams);


        gopts0 = zeros(nconds,ngams);


        rcond_max = 800;
        rcond_min = 1;


        for i=1:ngams
%
        prin2('gamma=',gams(i),1);

        prini(0,0);
        [gopts0(:,i),rconds] = errbetas_rconds0(gams(i),rcond_min,rcond_max,nconds);
        prini(13,0);
    end


        prin2('gopts0=',gopts0(:,1),nconds);
        prin2('gopts0=',gopts0(:,2),nconds);


%
%        plot optimal betas for each gamma
%
        hold on;
        box on;
        grid on;

        for i=1:ngams
%
        plot(rconds,gopts0(:,i),'linewidth',2);

        sss = int2str(log10(gams(i)));
        legs(i) = strcat("$\gamma=10^{",sss,"}$");

    end

        legend(legs,'interpreter','latex','fontsize',20,'location','southeast');
        set(figure(1),'Position',[500,500,900,800]);

        xlim([rcond_min-10,rcond_max]);
        ylim([min(gopts0(:,ngams)),max(gopts0(:,1))+.01]);


        xlabel('$\kappa$','interpreter','latex','fontsize',20);
        ylabel('$\beta$','interpreter','latex','fontsize',20);



        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'betas_errs_conds','-depsc','-r0')



        end        
%
%
%
%
%
        function [gopts0,rconds] = errbetas_rconds0(gam,rcond_min,rcond_max,nconds)

%
%        evaluates maximum beta where new method beats OptShrink,
%        for equispaced variances with specified condition numbers
%

        rconds = zeros(nconds,1);

        bbps = zeros(nconds,1);
        gopts0 = zeros(nconds,1);


        m=2000;


        prinf('m=',m,1);
        prin2('gam=',gam,1);


        for i=1:nconds
%
        prinf('i=',i,1);
        rconds(i) = rcond_min + (rcond_max-rcond_min)*(i-1)/(nconds-1);
    end


        prin2('conds=',rconds,nconds);


        rcondmax = rconds(nconds);
        vars_max = make_vars(rcondmax,m);

        prin2('rcondmax=',rcondmax,1);


        prinf('m=',m,1);




        for i=1:nconds
%
        vars = make_vars(rconds(i),m);
        xcut_bbp = eval_bbp(vars,m,gam);

        prin2('bbp=',xcut_bbp,1);


        sx = 1.1*xcut_bbp;



        prin2('vars=',vars,10);

%%%        prinstop;

        [gopts0(i),err_opt,dmu,tau,ierr] = maxbeta_err_secant(sx,vars,m,gam);

    end

%%%        prini(13,0);

        prin2('gopts0=',gopts0,nconds);



%%%        prinstop;
        end
%
%
%
%
%
        function test_maxbeta555()


        alp=12;
        m=1000;
        vars = make_vars2(alp,m);

        gam=.5;
        prinf('m=',m,1);

        nn = 1000;


        xcut_bbp = eval_bbp(vars,m,gam);

        prin2('bbp=',xcut_bbp,1);


        sx = 1.1*xcut_bbp;

        prin2('sx=',sx,1);


%%%        prinstop;

        [gopt0_dumb,gams0,errs,vals] = maxbeta_err_dumb(sx,vars,m,gam,nn);
%%%        prin2('err=',err,1);

        prin2('gopt0_dumb=',gopt0_dumb,1);



        [gopt0,err_opt,dmu,tau,ierr] = maxbeta_err_secant(sx,vars,m,gam);


        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gopt0,tau);
        val = err - err_opt;

        prin2('val=',val,1);

        prinf('ierr=',ierr,1);
        prinstop;


%
%        plot differences and root
%
        dzero=0;

        hold on;
        plot(gams0,vals);
        plot(gams0,dzero*ones(nn,1),'k');

        plot(gopt0,dzero,'r*');

        prinstop;


        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
            fshr_draw(sx,vars,m,n,k,n0,dds);


        dmu_hat = fshr_estmu(cov0,m);
        prin2('dmu_hat=',dmu_hat,1);

        prin2('xmat=',xmat,10);
        prin2('yf=',yf,10);

        chk0 = norm(xf+epf - yf);
        prin2('chk0=',chk0,1);


        prin2('dmu_hat=',dmu_hat,1);
        prin2('dmu=',dmu,1);




        prinstop;

%%%        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
%%%            fshr_pop2pars(sx,dmu,gam,gam0,tau);




%%%        prin2('err=',err,1);



        prinstop;


        iorth=0;
        [xhat,xfhat,errhat,khat,uyfw,vy,tvals,etas,taus_hat,sx_hat,couts,...
            cinns,dmu] = fshr_approx(y,cov0,gam0,m,n,k,iorth);

        errhat = errhat^2;


        err_true = norm(xhat-xmat,'fro')^2;

        prin2('err_true=',err_true,1);
        prin2('err=',err,1);
        prin2('errhat=',errhat,1);



        [err_opt,sy88,cout88,cinn88] = err_optshrink(sx,vars,gam,m);
        err_opt = err_opt^2;

        [xhat3,khat3,errhat3] = optshrink44(y,m,n,k,vars,ep);
        erropt_true = norm(xhat3-xmat,'fro')^2;


        prin2('err_opt=',err_opt,1);
        prin2('erropt_true=',erropt_true,1);


        prin2('xhat3=',xhat3,10);


        prinstop;
        end
%
%
%
%
%
        function [gopt0,err_opt,dmu,tau,ierr] = maxbeta_err_secant(sx,vars,m,gam)

        gams0 = zeros(1000,1);
        vals = zeros(1000,1);
        errs = zeros(1000,1);

        ierr=0;


        dmu = sum(vars) / m;

        dds=ones(m,1);
        k=1;
        tau = fshr_taus_exact(dds,vars,m,k);

        prin2('dmu=',dmu,1);
        prin2('tau=',tau,1);

%%%        prinstop;


        [err_opt,sy88,cout88,cinn88] = err_optshrink(sx,vars,gam,m);
        err_opt = err_opt^2;


        prin2('err_opt=',err_opt,1);



        ggg=1-sqrt(eps);
        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,ggg,tau);
        val = err - err_opt;

        prin2('val=',val,1);

        if (val <= 0 && err>0)
%
        gopt0=1;

        prin2('returning with gopt0=',gopt0,1);

        return;
    end


%%%        prinstop;


%
%        initial value, with positive error and value
%

        gleft0 = 0;
        gright0 = 1;

        gg0 = .25;

        for ijk=1:100

        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gg0,tau);
        val = err - err_opt;

        prin2('val=',val,1);
        prin2('err=',err,1);


%
%        if too far right, move left
%
        if (err == 0)
%
        gright0 = gg0;
        gg0 = (gg0 + gleft0) / 2;

        prinf('left on ijk=',ijk,1);
        prin2('gg0=',gg0,1);

        continue;
    end

%
%        if too far left, move right
%
        if (val<0 && err>0)
%
        gleft0 = gg0;
        gg0 = (gg0 + gright0) / 2;

        prinf('right on ijk=',ijk,1);
        prin2('gg0=',gg0,1);
        continue;
    end

%
%       check if done
%
        if (val >=0 && err>0)
        break;
    end
    end


        prin2('gg0=',gg0,1);

        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gg0,tau);
        val = err - err_opt;

        prin2('val=',val,1);
        prin2('err=',err,1);


%
%        find value between root and gg0
%

        hh0 = 0;
        for ijk=1:100
%
        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,hh0,tau);
        val = err - err_opt;

        if (val < 0)
%
        hh0 = (hh0 + gg0) / 2;
    end

    end

        prin2('hh0=',hh0,1);
        prin2('gg0=',gg0,1);


        [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,hh0,tau);
        val = err - err_opt;

        prin2('val=',val,1);

%%%        prinstop;

        gams0(1) = gg0;
        gams0(2) = hh0;


        [tval,eta,errs(1),syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gams0(1),tau);
        vals(1) = errs(1) - err_opt;

        [tval,eta,errs(2),syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gams0(2),tau);
        vals(2) = errs(2) - err_opt;

        prin2('vals=',vals,10);


        kplus=0;
        thresh = sqrt(eps)*100;

%
%        find root using secant method
%
        for i=3:100
%
        grad = (vals(i-1) - vals(i-2)) / (gams0(i-1) - gams0(i-2));

        gams0(i) = gams0(i-1) - vals(i-1)/grad;

%%%        xcuts(i) = fshr_xcut_tau(gam,gams0(i),tau);
%%%        vals(i) = xcuts(i) - xcut_bbp;

        [tval,eta,errs(i),syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gams0(i),tau);
        vals(i) = errs(i) - err_opt;

        if (abs(vals(i)) < thresh)
%
        kplus = kplus + 1;
    end

        if (kplus == 2)
%
        break;
    end

    end

        gopt0 = gams0(i);

        prin2('vals=',vals,10);
        prin2('gams0=',gams0,10);

        prin2('gopt0=',gopt0,1);


        if (kplus<2)

        prinf('failure to converge, kplus=',kplus,1);
        ierr=1;

    end

        if (errs(i)==0)

        prin2('bad region, value=',vals(i),1);
        prin2('bad region, err=',errs(i),1);
        ierr=ierr+2;
    end

%%%        prinstop;
        end
%
%
%
%
%
        function [gopt0,gams0,errs,vals,nn] = maxbeta_err_dumb(sx,vars,m,gam,nn)

%
%        evaluate errors as function of beta (gamm0),
%        and approximates maximum beta where errors cross
%
        [err_opt,sy88,cout88,cinn88] = err_optshrink(sx,vars,gam,m);
        err_opt = err_opt^2;


        prinf('nn=',nn,1);

        gams0 = zeros(nn,1);
        errs = zeros(nn,1);
        vals = zeros(nn,1);

        dmu = sum(vars) / m;

        dds=ones(m,1);
        k=1;
        tau = fshr_taus_exact(dds,vars,m,k);

%
%        evaluate errors on grid
%
        gmin0=0;
        gmax0=.999;

        for i=1:nn
%
        prinf('i=',i,1);

        gams0(i) = gmin0 + (gmax0 - gmin0)*(i-1)/(nn-1);

        [tval,eta,errs(i),syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gams0(i),tau);

        vals(i) = errs(i) - err_opt;

    end

        gopt0=1;
%
%        find approximate root
%
        for i=1:nn-1
%
        if (vals(i) <= 0 && vals(i+1) >= 0)
%
        gopt0 = (gams0(i) + gams0(i+1))/2;
        break;

    end

    end
        prin2('gams0=',gams0,nn);
        prin2('errs=',errs,nn);


        end
%
%
%
%
%
        function [m,n,n0,k,rcond,sx,vars,dds,gam,gam0,taus,xcut,ycut] = ...
            fshr_genpars555()
%
%        generates example problem parameters
%
        rcond = 1200;

        prin2('rcond=',rcond,1);

        m=1200;
        n=floor(2*m);
        n0=floor(3*m);

        k=1;
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

        vars = vars / rcond;
        prin2('vars=',vars,10);
        prin2('vars(m)=',vars(m),1);

%%%        prinstop;


%
%        signal vector weights
%
        dds=ones(m,k);

        for j=1:k
%
        dnorm = norm(dds(:,j));
        dds(:,j) = dds(:,j) / dnorm * sqrt(m);
    end

        prin2('dd1=',dds(:,1),10);

        chk0 = abs(norm(dds(:,1)) - sqrt(m));
        prin2('chk0=',chk0,1);


        taus = fshr_taus_exact(dds,vars,m,k);
        vars = vars*taus;


        taus = fshr_taus_exact(dds,vars,m,k);


        prin2('taus=',taus,k);


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
        sx(i) = xcut*(1.12 + (k-i+1)/k) / sqrt(taus(i));


    end

        prin2('sx=',sx,k);
        prin2('threshold=',xcut / sqrt(taus(1)),1);


        randn(1,101);
        randn(1,101);


        end
%
%
%
%
%
        function plot_thresh()
%
%        minimum detectable signals as function of beta/gamm0,
%        for different gamma values
%
        n=5000;

        ngams=8;

        gams0 = zeros(n,1);
        gams = zeros(ngams,1);
        xcuts = zeros(n,ngams);


        gmin = 0;
        gmax = .6;


        pp = [3,2,1,0,-1,-3,-5];


        for i=1:ngams-1
%
        gams(i) = 2^pp(i);
    end

        gams(ngams) = 0;

        prin2('gams=',gams,ngams);

        for i=1:n
%
        gams0(i) = gmin*(n-i)/(n-1) + gmax*(i-1)/(n-1);
%

        for j=1:ngams
%
        [ycut,xcuts(i,j)] = fshr_cuts(gams(j),gams0(i));
    end
    end

        prin2('gam=',gams(1),1);
        prin2('gams0=',gams0,n);


        prin2('xcuts=',xcuts,n);


        dzero=0;
        [ycut,xcut0] = fshr_cuts(gams(1),dzero);

        prin2('xcut0=',xcut0,1);
%
        figure(1);
        hold on;
        box on;
        grid on;
%
%%%        plot(gams0,xcuts(1:n,4),'k','Linewidth',2);
%%%        plot(gams0,xcuts(1:n,2),'r--','Linewidth',2);
%%%        plot(gams0,xcuts(1:n,3),'g-.','Linewidth',2);
%%%        plot(gams0,xcut0*ones(n,1),'Color','r','Linewidth',2);


        for j=1:ngams

        plot(gams0,xcuts(:,j),'Linewidth',2);
    end

        xlabel('$\beta$','interpreter','latex','fontsize',20);
        ylabel('$\sigma_{\mathrm{thresh}}$','interpreter','latex','fontsize',20);

        ymin = 0;
        ymax = xcuts(n,1);

        ylim([ymin,ymax]);

%%%        prinstop;
%%%        yts=linspace(ymin,ymax,5);
%%%        yticks(yts);


        for j=1:ngams-1
%
        ll = "$\gamma = 2^{"
        ll = strcat(ll,int2str(pp(j)))
        ll = strcat(ll,"}$")

        legs(j) = ll
    end

        legs(ngams) = "$\gamma=0$"

%%%        legend({'$\gamma=3$','$\gamma=1$','$\gamma=1/10000$'},...
%%%            'interpreter','latex','fontsize',18,'location','northwest')

%%%        legend(legs,'interpreter','latex','location','northwest','fontsize',20)
        legend(legs,'interpreter','latex','location','southeast','fontsize',20)


%%%        set(figure(1),'position',[500,500,1000,870]);

        set(figure(1),'Position',[500,500,900,800]);


        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'thresholds','-depsc','-r0')



        prinstop;
        end
%
%
%
%
%
        function ttt()

        randn(100,1);
        randn(100,1);

        nplot=200;

        m=1000;

        k=1;

        rconds = [1,100,1000];
        gams = [2,1,0.5];


        errs_mat = zeros(nplot,9);
        pars = zeros(9,2);

        errs2 = zeros(nplot,1);
        errs3 = zeros(nplot,1);


        prin2('gams=',gams,3);
%%%        prinstop;


%
%        detectable signal singular value
%
        gmin0 = 0;
        gmax0 = .95;



        gmax0 = .5;
        gmax0 = .3;


        rcond_min = min(rconds);

        rcond_min = 15;


        gam_max = max(gams);

        gam_max=1;

        vars = make_vars(rcond_min,m);

        sig=2;
        vars=sig^2*vars;

        prin2('vars=',vars,100);


        [xcut,xcut2,xcut3] = cutoff_both(vars,m,gam_max,gmax0);

        sx = 1.5*xcut3;



        ccc = sig*gam_max^(1/4);

        prin2('ccc=',ccc,1);


        prin2('xcut2=',xcut2,1);
        prin2('xcut3=',xcut3,1);



        prinstop;


        k=3;
        dds=ones(m,k);
        n=floor(m/gam_max);
        n0=floor(m/gmax0);

        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
            fshr_draw(sx,vars,m,n,k,n0,dds);

        ifig=1;
        fshr_wacompare(yf,gam_max,gmax0,ifig);

        ifig=2;
        fshr_wacompare(y,gam_max,gmax0,ifig);


        prinstop;



        prin2('sx=',sx,1);

%
%        signal singular vector weights
%
        dds=ones(m,1);
%
        prin2('dds=',dds,10);


        gams0 = zeros(nplot,1);

        for i=1:nplot
%
%%%        prinf('notch i=',i,1);

        gams0(i) = gmin0 + (gmax0-gmin0)*(i-1)/(nplot-1);
    end




        ijk=1;

        for j=1:3

        for i=1:3

        rcond = rconds(i);
%
%        noise variances and tau
%
        vars = make_vars(rcond,m);
%
        tau = fshr_taus_exact(dds,vars,m,k);
        prin2('tau=',tau,1);


        gam = gams(j);

%%%        prin2('vars=',vars,m);
%%%        prin2('vars end=',vars(m),1);


        [errs,errs2(ijk),errs3(ijk),gams0] = errs_all(sx,vars,gam,gmin0,gmax0,m,nplot);

        errs_mat(:,ijk) = errs;
        pars(ijk,1) = rconds(i);
        pars(ijk,2) = gams(j);

        ijk=ijk+1;
    end

    end
        prin2('errs=',errs,nplot);
        prin2('err2=',errs2,9);
        prin2('err3=',errs3,9);



        prin2('parameter, condition numbers=',pars(:,1),9)

        prin2('parameter, gammas=',pars(:,2),9)



%%%        hold on;
%%%        ijk=1;
%%%        plot(gams0,errs_mat(:,ijk));
%%%        plot(gams0,errs2(ijk)*ones(nplot,1));
%%%        plot(gams0,errs3(ijk)*ones(nplot,1));

%%%        prinstop;
%
%        make plots
%
        figure(1);

        for ijk=1:9

        subplot(3,3,ijk);

        rcond = pars(ijk,1);
        gam = pars(ijk,2);

        prin2('rcond=',rcond,1);
        prin2('gam=',gam,1);




        hold on;
        box on;
        grid on;

        plot(gams0,errs_mat(:,ijk),'k','LineWidth',1);
        plot(gams0,errs3(ijk)*ones(nplot,1),'r-.','LineWidth',1);
        plot(gams0,errs2(ijk)*ones(nplot,1),'b--','LineWidth',1);

        xlim([min(gams0),max(gams0)]);
%
        ymin=.9*min(errs2(ijk),errs3(ijk));
        ymax=1.01*max([max(errs_mat(:,ijk)),errs2(ijk),errs3(ijk)]);
        ylim([ymin,ymax]);


%%%        prinstop;

%
%        make plot title
%
        if (gam < 1)
%
        tstr1 = strcat(['$\gamma=$',' 1/',int2str(1/gam),',']);
    end
        if (gam >= 1)
%
        tstr1 = strcat(['$\gamma=$',' ',int2str(gam),',']);
    end

        tstr2 = strcat([' $\kappa= $ ',' ',int2str(rcond)]);
        tstr=strcat(tstr1,tstr2);
        title(tstr,'interpreter','latex','fontsize',12);

        xlabel('$\beta$','interpreter','latex','fontsize',12);
        ylabel('Estimated AMSE','interpreter','latex','fontsize',12);

    end


        set(figure(1),'position',[500,500,1100,900]);


%
%        add legend, resize, and save
%
        subplot(3,3,1);
        legend({'Pseudo-whitening','No whitening','Oracle whitening'},...
            'interpreter','latex','fontsize',12,'location','northwest');

%
%        save figure to file
%
        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'plot_errs','-depsc','-r0')





        prinstop;



        n=floor(m/gam);
        n0=floor(m/gam0);

        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
            fshr_draw(sx,vars,m,n,k,n0,dds);

        [err_real,err2_real,err3_real,err4_real,err5_real,xhat,errhat,khat,xhat2,...
            khat2,xhat3,khat3,errhat3,xhat4,khat4,xhat5,khat5] = ...
            allmethods(y,xmat,cov0,ep0,gam0,vars,ep,m,n,k);

        prin2('err3, real=',err3_real,1);
        prin2('err3=',err3,1);



        prin2('err4, real=',err4_real,1);
        prin2('err2=',err2,1);



        end
%
%
%
%
%
        function plot_errors()

        randn(100,1);
        randn(100,1);

        nplot=200;

        m=1000;

        k=1;

        rconds = [1,5000,50000];
        gams = [2,1,0.5];


        errs_mat = zeros(nplot,9);
        pars = zeros(9,2);

        errs2 = zeros(nplot,1);
        errs3 = zeros(nplot,1);


        prin2('gams=',gams,3);
%%%        prinstop;


%
%        detectable signal singular value
%
        gmin0 = 0;
        gmax0 = .95;

        rcond_min = min(rconds);
        rcond_max = max(rconds);
        gam_max = max(gams);
        vars = make_vars(rcond_max,m);
        xcut = cutoff_both(vars,m,gam_max,gmax0);

        sx = xcut+2;



        prin2('sx=',sx,1);

%
%        signal singular vector weights
%
        dds=ones(m,1);
%
        prin2('dds=',dds,10);


        gams0 = zeros(nplot,1);

        for i=1:nplot
%
%%%        prinf('notch i=',i,1);

        gams0(i) = gmin0 + (gmax0-gmin0)*(i-1)/(nplot-1);
    end




        ijk=1;

        for j=1:3

        for i=1:3

        rcond = rconds(i);
%
%        noise variances and tau
%
        vars = make_vars(rcond,m);
%
        tau = fshr_taus_exact(dds,vars,m,k);
        prin2('tau=',tau,1);


        gam = gams(j);

%%%        prin2('vars=',vars,m);
%%%        prin2('vars end=',vars(m),1);


        [errs,errs2(ijk),errs3(ijk),gams0] = errs_all(sx,vars,gam,gmin0,gmax0,m,nplot);

        errs_mat(:,ijk) = errs;
        pars(ijk,1) = rconds(i);
        pars(ijk,2) = gams(j);

        ijk=ijk+1;
    end

    end
        prin2('errs=',errs,nplot);
        prin2('err2=',errs2,9);
        prin2('err3=',errs3,9);



        prin2('parameter, condition numbers=',pars(:,1),9)

        prin2('parameter, gammas=',pars(:,2),9)



%%%        hold on;
%%%        ijk=1;
%%%        plot(gams0,errs_mat(:,ijk));
%%%        plot(gams0,errs2(ijk)*ones(nplot,1));
%%%        plot(gams0,errs3(ijk)*ones(nplot,1));

%%%        prinstop;
%
%        make plots
%
        figure(1);

        for ijk=1:9

        subplot(3,3,ijk);

        rcond = pars(ijk,1);
        gam = pars(ijk,2);

        prin2('rcond=',rcond,1);
        prin2('gam=',gam,1);


        hold on;
        box on;
        grid on;

        plot(gams0,errs_mat(:,ijk),'k','LineWidth',1);
        plot(gams0,errs3(ijk)*ones(nplot,1),'r-.','LineWidth',1);
        plot(gams0,errs2(ijk)*ones(nplot,1),'b--','LineWidth',1);

        xlim([min(gams0),max(gams0)]);
%
        ymin=.9*min(errs2(ijk),errs3(ijk));
        ymax=1.05*max([max(errs_mat(:,ijk)),errs2(ijk),errs3(ijk)]);
        ylim([ymin,ymax]);


%%%        prinstop;

%
%        make plot title
%
        if (gam < 1)
%
        tstr1 = strcat(['$\gamma=$',' 1/',int2str(1/gam),',']);
    end
        if (gam >= 1)
%
        tstr1 = strcat(['$\gamma=$',' ',int2str(gam),',']);
    end

        tstr2 = strcat([' $\kappa= $ ',' ',int2str(rcond)]);
        tstr=strcat(tstr1,tstr2);
        title(tstr,'interpreter','latex','fontsize',12);

        xlabel('$\beta$','interpreter','latex','fontsize',12);
        ylabel('Estimated AMSE','interpreter','latex','fontsize',12);

    end


        set(figure(1),'position',[500,500,1100,900]);


%
%        add legend, resize, and save
%
        subplot(3,3,1);
        legend({'Pseudo-whitening','No whitening','Oracle whitening'},...
            'interpreter','latex','fontsize',12,'location','northwest');

%
%        save figure to file
%
        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'plot_errs','-depsc','-r0')





        prinstop;



        n=floor(m/gam);
        n0=floor(m/gam0);

        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
            fshr_draw(sx,vars,m,n,k,n0,dds);

        [err_real,err2_real,err3_real,err4_real,err5_real,xhat,errhat,khat,xhat2,...
            khat2,xhat3,khat3,errhat3,xhat4,khat4,xhat5,khat5] = ...
            allmethods(y,xmat,cov0,ep0,gam0,vars,ep,m,n,k);

        prin2('err3, real=',err3_real,1);
        prin2('err3=',err3,1);



        prin2('err4, real=',err4_real,1);
        prin2('err2=',err2,1);



        end
%
%
%
%
%
        function plot_errors2()

        randn(100,1);
        randn(100,1);

        nplot=200;

        m=1000;

        k=1;

        rconds = [1,100,1000];


%%%        alps = [1,4,6];
        alps = [1,3,5];

        gams = [2,1,0.5];


        errs_mat = zeros(nplot,9);
        pars = zeros(9,2);

        errs2 = zeros(nplot,1);
        errs3 = zeros(nplot,1);


        prin2('gams=',gams,3);
%%%        prinstop;


%
%        detectable signal singular value
%
        gmin0 = 0;
        gmax0 = .95;


%%%        alp0 = min(alps);
        alp0 = max(alps);
        vars = make_vars2(alp0,m);


%%%        prin2('vars=',vars,10);
%%%        prinstop;

        gam_max = max(gams);
        [xcut,xcut2,xcut3] = cutoff_both(vars,m,gam_max,gmax0);

        sx = xcut+1;

        prin2('gam_max=',gam_max,1);
        prin2('alp0=',alp0,1);

%%%        prinstop;


        prin2('xcut2=',xcut2,1);
        prin2('xcut3=',xcut3,1);
        prin2('xcut=',xcut,1);

        prin2('sx=',sx,1);

%%%        prinstop;

%
%        signal singular vector weights
%
        dds=ones(m,1);
%
        prin2('dds=',dds,10);


        gams0 = zeros(nplot,1);

        for i=1:nplot
%
%%%        prinf('notch i=',i,1);

        gams0(i) = gmin0 + (gmax0-gmin0)*(i-1)/(nplot-1);
    end


        ijk=1;

        for j=1:3

        for i=1:3

        alp = alps(i);
%
%        noise variances and tau
%
        vars = make_vars2(alp,m);
%
        tau = fshr_taus_exact(dds,vars,m,k);
        prin2('tau=',tau,1);


        gam = gams(j);

%%%        prin2('vars=',vars,m);
%%%        prin2('vars end=',vars(m),1);


        [errs,errs2(ijk),errs3(ijk),gams0] = errs_all(sx,vars,gam,gmin0,gmax0,m,nplot);

        errs_mat(:,ijk) = errs;
        pars(ijk,1) = alps(i);
        pars(ijk,2) = gams(j);

        ijk=ijk+1;
    end

    end
        prin2('errs=',errs,nplot);
        prin2('err2=',errs2,9);
        prin2('err3=',errs3,9);



        prin2('parameter, alphas=',pars(:,1),9)

        prin2('parameter, gammas=',pars(:,2),9)



%%%        hold on;
%%%        ijk=1;
%%%        plot(gams0,errs_mat(:,ijk));
%%%        plot(gams0,errs2(ijk)*ones(nplot,1));
%%%        plot(gams0,errs3(ijk)*ones(nplot,1));

%%%        prinstop;
%
%        make plots
%
        figure(1);

        for ijk=1:9

        subplot(3,3,ijk);

        alp = pars(ijk,1);
        gam = pars(ijk,2);

        prin2('rcond=',alp,1);
        prin2('gam=',gam,1);




        hold on;
        box on;
        grid on;

        plot(gams0,errs_mat(:,ijk),'k','LineWidth',1);
        plot(gams0,errs3(ijk)*ones(nplot,1),'r-.','LineWidth',1);
        plot(gams0,errs2(ijk)*ones(nplot,1),'b--','LineWidth',1);

        xlim([min(gams0),max(gams0)]);
%
        ymin=.85*min(errs2(ijk),errs3(ijk));
%%%        ymin=0;
        ymax=1.04*max([max(errs_mat(:,ijk)),errs2(ijk),errs3(ijk)]);

        ymin
        ymax

        errs3

        ylim([ymin,ymax]);



%%%        prinstop;

%
%        make plot title
%
        if (gam < 1)
%
        tstr1 = strcat(['$\gamma=$',' 1/',int2str(1/gam),',']);
    end
        if (gam >= 1)
%
        tstr1 = strcat(['$\gamma=$',' ',int2str(gam),',']);
    end

        tstr2 = strcat([' $\alpha= $ ',' ',int2str(alp)]);
        tstr=strcat(tstr1,tstr2);
        title(tstr,'interpreter','latex','fontsize',12);

        xlabel('$\beta$','interpreter','latex','fontsize',12);
        ylabel('Estimated AMSE','interpreter','latex','fontsize',12);

    end


        set(figure(1),'position',[500,500,1100,900]);


%
%        add legend, resize, and save
%
        subplot(3,3,1);
        legend({'Pseudo-whitening','No whitening','Oracle whitening'},...
            'interpreter','latex','fontsize',12,'location','northwest');

%
%        save figure to file
%
        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'plot_errs2','-depsc','-r0')





        prinstop;



        n=floor(m/gam);
        n0=floor(m/gam0);

        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
            fshr_draw(sx,vars,m,n,k,n0,dds);

        [err_real,err2_real,err3_real,err4_real,err5_real,xhat,errhat,khat,xhat2,...
            khat2,xhat3,khat3,errhat3,xhat4,khat4,xhat5,khat5] = ...
            allmethods(y,xmat,cov0,ep0,gam0,vars,ep,m,n,k);

        prin2('err3, real=',err3_real,1);
        prin2('err3=',err3,1);



        prin2('err4, real=',err4_real,1);
        prin2('err2=',err2,1);



        end
%
%
%
%
%
        function vars = make_vars2(alp,m)

        vars=zeros(m,1);
        tts=zeros(m,1);

        for i=1:m
%
        tts(i) = 2*(i-1)/(m-1) + 1;
        vars(i) = tts(i)^(-alp);
%%%        vars(i) = tts(i)^(alp);

    end

%%%        return;

%
%        rescale variances so that tau is 1
%
        dds = ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);

        vars = vars * tau;

%%%        return;


        prin2('tau=',tau,1);
%%%        prin2('vars=',vars,m);

%%%        prinstop;

%
%        check that tau is 1
%
        tau = fshr_taus_exact(dds,vars,m,1);
        prin2('tau=',tau,1);

        chk0 = tau-1;
        prin2('chk0=',chk0,1);

%%%        prinstop;


        return;

        prin2('tts=',tts,m);
        prin2('vars=',vars,m);

%%%        plot(tts,vars);

%%%        prin2('ss3=',sqrt(3),1);

%%%        prinstop;

%%%        vars = vars / sum(vars);

%%%        vars = vars / alp;





 
        end
%
%
%
%
%
        function [xcut,xcut2,xcut3] = cutoff_both(vars,m,gam,gam0)

        dds=ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);
        xcut2 = fshr_xcut_tau(gam,gam0,tau);

        n=1;
        dzero = 0;
        awhts=ones(1,m) / m;
        bwhts=1;
        bs=1;

        prinf('n=',n,1);
        prinf('m=',m,1);
        prin2('gam=',gam,1);
%%%        prinstop;
        [ell,bedge,err] = mpbdry_thresh(vars,bs,awhts,bwhts,m,n,gam);

        xcut3 = sqrt(ell);

        xcut = max(xcut2,xcut3);


        prin2('xcut=',xcut,1);


        prin2('xcut2=',xcut2,1);
        prin2('xcut3=',xcut3,1);
        prin2('gam0=',gam0,1);

%%%        prinstop;


        end
%
%
%
%
%
        function [errs,err2,err3,gams0] = errs_all(sx,vars,gam,gmin0,gmax0,m,nplot)

        errs = zeros(nplot,1);
%
%        parameter mu
%
        dmu=0;
        for i=1:m
%
        dmu = dmu+vars(i)/m;
    end

        dds=ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);

        prin2('dmu=',dmu,1);
        prin2('tau=',tau,1);


        err2 = err_oracle(sx,tau,dmu,gam);
        prin2('err2=',err2,1);



        [err3,sy88,cout88,cinn88] = err_optshrink(sx,vars,gam,m);


        prin2('err3=',err3,1);


        errs = zeros(nplot,1);


        gams0 = zeros(nplot,1);

        for i=1:nplot
%
%%%        prinf('notch i=',i,1);

        gams0(i) = gmin0 + (gmax0-gmin0)*(i-1)/(nplot-1);


        [tval,eta,errs(i),syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gams0(i),tau);
        errs(i) = sqrt(errs(i));

        prin2('error=',errs(i),1);

    end

        prin2('gams0=',gams0,nplot);
        prin2('errs=',errs,nplot);

%%%        prinstop;

        prin2('err2=',err2,1);
        prin2('err3=',err3,1);

        prin2('gams0=',gams0,nplot);

        prin2('sx=',sx,1);

%%%        prinstop;

        end
%
%
%
%
%
        function [err,sy,cout,cinn] = err_optshrink(sx,vars,gam,m)
%
%%%        prin2('sx=',sx,1);
%%%        prin2('vars=',vars,10);

%%%        gam=m/n;
        ell=sx^2;

        n=1;

        bs=ones(1,n);
        awhts=ones(1,m) / m;
        bwhts=ones(1,n) / n;
        [rlam,cout,cinn] = mpbdry_sforw(ell,vars,bs,awhts,bwhts,m,n,gam);

        sy=sqrt(rlam);

        err=ell*(1-cout^2*cinn^2);
        err=sqrt(err);


        return;

        prin2('err=',err,1);

        prin2('rlam=',rlam,1);
        prin2('cout=',cout,1);
        prin2('cinn=',cinn,1);

        return;

        bedge = mpbdry_edge(vars,bs,awhts,bwhts,m,n,gam);
        prin2('ell=',ell,1);
        prin2('bedge=',bedge,1);
        prinstop;

        end
%
%
%
%
%
        function plot_shrinkers()
%
        randn(1,1000);
        randn(1,1000);
        randn(1,1000);
        randn(1,1000);
        randn(1,1000);
        randn(1,1000);
        randn(1,1000);

        m=1000;
        n=floor(2*m);

        k=1;

        rcond = 500;
        prinf('rcond=',rcond,1);
%
        vars = zeros(m,1);

        gam=m/n;
        prin2('gam=',gam,1);

%
%        noise variances
%
        for i=1:m
%
        vars(i) = rcond*(i-1)/(m-1) + (m-i)/(m-1);
    end
%
%%%        vars = linspace(1,rcond,m)';

        vars = vars/rcond;
        prin2('vars=',vars,10);
        prin2('vars end=',vars(m),1);



%%%        prinstop;

%
%        signal singular vector weights
%
        dds=ones(m,1);
%
        prin2('dds=',dds,m);

%
%        exact values of tau
%
        tau = fshr_taus_exact(dds,vars,m,k);
        prin2('tau=',tau,1);

%%%        prinstop;

%
%        parameter mu
%
        dmu=0;
        for i=1:m
%
        dmu = dmu+vars(i)/m;
    end

        prin2('dmu=',dmu,1);

%%%        prinstop;


%%%        tic;

%%%        [ycut,xcut] = allcuts(gam,0);
        [ycut,xcut] = fshr_cuts(gam,0);
        xcut=xcut/sqrt(tau);
        smin = xcut-.1;
        smax = xcut + 5;

        prin2('smin=',smin,1);
        prin2('smax=',smax,1);


        nvals=500;

        ngams0 = 4;
        gams0 = zeros(ngams0,1);

        gams0(1) = 0;
        gams0(2) = 0.3;
        gams0(3) = 0.5;
        gams0(4) = 0.7;


        tvals = zeros(nvals,ngams0);
        errs = zeros(nvals,ngams0);

        xvals = zeros(nvals,1);
        yvals = zeros(nvals,ngams0);


        for i=1:nvals
%
        xvals(i) = smin + (smax-smin)*(i-1)/(nvals-1);


        for j=1:ngams0
%
        prini(0,0);
        [tvals(i,j),errs(i,j),eta,unorm] = pop2shrinker(xvals(i),tau,dmu,gam,gams0(j));
        prini(13,0);
    end

    end


%
%        for each gam0, map population values to empirical values
%
        for j=1:ngams0

        [ycut2,xcut2] = fshr_cuts(gam,gams0(j));
        xcut2=xcut2/sqrt(tau);

        for i=1:nvals
%
%%%        [yvals(i,j),cinn,cw,cc,vpsi,dpsi,sbar,sbder,ss,sder,res,der] = ...
%%%            pop2emp555(sqrt(tau)*xvals(i),gam,gams0(j),tau);

        if (xvals(i) <= xcut2)
%
        yvals(i,j) = ycut2;
        continue;
    end

        yvals(i,j) = fshr_spikforw(xvals(i),gam,gams0(j),tau);
    end

    end

        prin2('yvals=',yvals(:,j),10);

        prin2('xvals=',xvals,10);

        j=4;

        n=floor(m/gam);
        n0=floor(m/gams0(j));
        k=1;

        prinf('m=',m,1);
        prinf('n=',n,1);
        prinf('n0=',n0,1);


        if (0)
%
        i=500;
        sx=xvals(i);

        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
            fshr_draw(sx,vars,m,n,k,n0,dds);


        prin2('yvals(i,j)=',yvals(i,j),1);
        prin2('syf=',syf,1);

        qq = sx/sqrt(1 - gams0(j)) * sqrt(tau);
        prin2('qq=',qq,1);

        prinstop;
    end

        nvals2=floor(nvals/2);
        plot_x2t(xvals(1:nvals2),tvals(1:nvals2,:),nvals2,gams0,ngams0);


        plot_y2t(yvals,tvals,nvals,gams0,ngams0);




        prinstop;



        end
%
%
%
%
%
        function [tval,err,eta,unorm] = pop2shrinker(sx,tau,dmu,gam,gam0)
%
%        evaluates optimal shrinker and error from population parameters
%
        tval=0;
        err=0;
        eta=0;
        unorm=0;

        [ycut,xcut] = fshr_cuts(gam,gam0);

        if (sx <= xcut / sqrt(tau))
%
        prinf('in shrinker, tval=',tval,1);
        return;
    end

        [sxf,syf,cinn,cw,cc,vpsi,dpsi,sbar,sbder,ss,sder,res,der] ...
             = pop2fmlas(sx,gam,gam0,tau);

        yy=syf^2;

        [upone,dupone] = fshr_upone_fmla(yy,gam,gam0);
        uptwo = fshr_uptwo_fmla(yy,gam,gam0);

        eone = -gam*res*(upone + yy*dupone) / dpsi + yy*sbar*(uptwo-ss^2)/dpsi;
        etwo = yy*sbar*ss^2 / dpsi;

        unorm = sqrt(etwo/tau + eone*dmu);
%
%        singular value, generalized value, squared error
%
        tval = sx*cinn*cw / unorm;
        eta = sx*cinn*cw / unorm^2;
        err = sx^2*(1 - cinn^2*cw^2 / unorm^2);

        err=sqrt(err);

        end
%
%
%
%
%
        function [sxf,syf,cinn,cw,cc,vpsi,dpsi,sbar,sbder,ss,sder,res,der] ...
             = pop2fmlas(sx,gam,gam0,tau)
%
        xval=sqrt(tau)*sx;
        syf = fshr_spikforw0(xval,gam,gam0);
        sxf = xval/sqrt(1-gam0);
%
        [cinn,cw,cc,vpsi,dpsi,sbar,sbder,ss,sder,res,der] = ...
            fshr_fmlas(syf,gam,gam0,tau);

        end
%
%
%
%
%
        function [cinn,cw,cc,vpsi,dpsi,sbar,sbder,ss,sder,res,der] = ...
            fshr_fmlas(yval,gam,gam0,tau)
%
        yy=yval^2;

        [vpsi,dpsi] = fshr_evalpsi(yy,gam,gam0);
        [sbar,sbder] = fshr_wachtbarst(yy,gam,gam0);

        [ss,sder] = fshr_wachtstiel_right(yy,gam,gam0);
        [res,der] = fshr_resolv(yy,gam,gam0);

        cinn = sbar*vpsi / dpsi;
        cinn = sqrt(cinn);

        cc = tau*res*vpsi/dpsi;
        cc=sqrt(cc);

%%%        cw = yy * ss^2 * xx * cinn^2;
        cw = yy * sbar * ss^2 / tau / dpsi;
        cw = sqrt(cw);


        prin2('sbar, inside fmlas=',sbar,1);
        prin2('ss squared, inside fmlas=',ss^2,1);
        prin2('tau, inside fmlas=',tau,1);
        prin2('dpsi, inside fmlas=',dpsi,1);
        prin2('yy, inside fmlas=',yy,1);




        return;

        cc2 = sqrt(res/dpsi/xx);

        chk0 = cc2-cc;
        prin2('chk0=',chk0,1);
        prinstop;

        end
%
%
%
%
%
        function plot_y2t(yvals,tvals,nvals,gams0,ngams0)

        figure(2);

%
%
%        plot as function of population value
%
        hold on;
        box on;
        grid on;
%

%
%        legend labels
%
        legs(1) = "$\beta=0$ (oracle)";

        for j=2:ngams0
%
        tstr = gam0_str(gams0(j));
        legs(j) = tstr;
    end

        xlabel('Observed singular value','interpreter','latex','fontsize',18);
        ylabel('Shrunken singular value','interpreter','latex','fontsize',18);

%
%        marker styles
%
        marks=["k","g-.","b--","r:"];

        for j=1:ngams0
%
        plot(yvals(:,j),tvals(:,j),marks(j),'LineWidth',2);
    end


%
%        black out the x axis
%
%%%        plot(yvals,zeros(1,nvals),'k');

        xlim([min(yvals(:)),max(yvals(:))]);
        ylim([0,tvals(nvals,1)]);


        legend(legs,'interpreter','latex','fontsize',18,'location','southeast')

        set(figure(2),'position',[500,500,700,600]);


%%%        prinstop;

        set(figure(2),'PaperPositionMode','auto')
        print(figure(2),'plot_y2t','-depsc','-r0')


%%%        prinstop;



        end
%
%
%
%
%
        function plot_x2t(xvals,tvals,nvals,gams0,ngams0)

        figure(1);

%
%
%        plot as function of population value
%
        hold on;
        box on;
        grid on;
%

%
%        legend labels
%
        legs(1) = "$\beta=0$ (oracle)";

        for j=2:ngams0
%
        tstr = gam0_str(gams0(j));
        legs(j) = tstr;
    end


        xlabel('Population singular value','interpreter','latex','fontsize',18);
        ylabel('Shrunken singular value','interpreter','latex','fontsize',18);

%
%        marker styles
%
        marks=["k","g-.","b--","r:"];

        for j=1:ngams0
%

        iii = tvals(:,j) > 0;

        plot(xvals(iii),tvals(iii,j),marks(j),'LineWidth',2);
    end


        xlim([xvals(1),xvals(nvals)]);
        ylim([0,tvals(nvals,1)]);


        legend(legs,'interpreter','latex','fontsize',18,'location','northwest')

        set(figure(1),'position',[500,500,700,600]);


%%%        prinstop;

        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'plot_x2t','-depsc','-r0')






        end
%
%
%
%
%
        function tstr = gam0_str(gam0)
%
        tstr = strcat(['$\beta=$',' 0.',int2str(10*gam0)]);

        end
%
%
%
%
%
        function plot_condn()

        load('compare_cond.mat');
        who;


        figure(1);
        hold on;
%%%        plot(log(rconds),rerrs_mean,'*-');
%%%        plot(log(rconds),rerrs2_mean,'*-');
%%%        plot(log(rconds),rerrs3_mean,'*-');
%%%        plot(log(rconds),rerrs4_mean,'*-');

        grid on;
        box on;

        lw=2;
        msize=10;

        plot(log2(rconds),log2(rerrs2_mean),'g>-','Linewidth',lw,'MarkerSize',msize);
%%%        plot(log2(rconds),log2(rerrs5_mean),'m^-','Linewidth',lw,'MarkerSize',msize);
        plot(log2(rconds),log2(rerrs_mean),'ks-','Linewidth',lw,'MarkerSize',msize);
        plot(log2(rconds),log2(rerrs3_mean),'ro-','Linewidth',lw,'MarkerSize',msize);
        plot(log2(rconds),log2(rerrs4_mean),'bd-','Linewidth',lw,'MarkerSize',msize);

%%%        legend('Plugin','Plugin, rank','Pseudo-whitening','OptShrink','Oracle whitening',...
%%%            'interpreter','latex','Location','Northwest','Fontsize',20)

        legend('Oracle with plugin','Pseudo-whitening','OptShrink','Oracle whitening',...
            'interpreter','latex','Fontsize',20)

        xlabel('$\log_{2}$ condition number','interpreter','latex','Fontsize',20);
        ylabel('$\log_{2}$ mean error','interpreter','latex','Fontsize',20);

        xlim([log2(rconds(1)),log2(rconds(nconds))]);

        ymin = min([log2(rerrs_mean(1)),log2(rerrs2_mean(1)),log2(rerrs3_mean(1)),...
            log2(rerrs4_mean(1)),]);
        ymax = max([log2(rerrs_mean(nconds)),log2(rerrs2_mean(nconds)),log2(rerrs3_mean(nconds)),...
            log2(rerrs4_mean(nconds)),]);
%%%        ylim([ymin,ymax]);

        set(figure(1),'Position',[500,500,900,800]);

%
%        save figure to file
%
        set(figure(1),'PaperPositionMode','auto')
        print(figure(1),'compare_known_rank','-depsc','-r0')



        prinstop;
        end
%
%
%
%
%
        function monte_condnums()
%
        iplot=0;

        if (iplot==1)
%
        plot_condn();
        prinstop;
    end


        nconds=12;
        ndraws=3;
%
        rerrs=zeros(nconds,ndraws);
        rerrs2=zeros(nconds,ndraws);
        rerrs3=zeros(nconds,ndraws);
        rerrs4=zeros(nconds,ndraws);
        rerrs5=zeros(nconds,ndraws);
%
        khats=zeros(nconds,ndraws);
        khats2=zeros(nconds,ndraws);
        khats3=zeros(nconds,ndraws);
        khats4=zeros(nconds,ndraws);
        khats5=zeros(nconds,ndraws);
%
        khat_maxes = zeros(nconds,5);
        khat_mins = zeros(nconds,5);
        khat_means = zeros(nconds,5);
        khat_medians = zeros(nconds,5);

%
%        set constant parameters
%
        m=600;
        n=floor(m);
        n0=floor(2*m);

        k=5;
        keep=5;

%%%        m=600;
%%%        n=floor(5*m);
%%%        n0=floor(2*m);
%%%
%%%        k=5;
%%%        keep=k;
%

        rconds = zeros(nconds,1);

%
%        noise condition numbers
%
        rconds(1) = 1;

        for i=2:nconds
%
        rconds(i) = 2*rconds(i-1);
    end

%
%        signal vector weights
%
        dds=ones(m,k);


        tic;

%
%        Monte Carlo draws
%
%%%        parfor iii=1:ndraws
        for iii=1:ndraws
%
        rng(iii);

        prinf('iii=',iii,1);

        [rerrs(:,iii),rerrs2(:,iii),rerrs3(:,iii),rerrs4(:,iii),rerrs5(:,iii)...
            khats(:,iii),khats2(:,iii),khats3(:,iii),khats4(:,iii),khats5(:,iii)] = ...
            sweep_condnums(rconds,dds,m,n,n0,k,keep,nconds);

    end


        time_mc = toc;
%%%        prin2('time, seconds=',time_mc,1);

        prinar2('rerrs=',rerrs,nconds,ndraws);
        prinar2('rerrs2=',rerrs2,nconds,ndraws);
        prinar2('rerrs3=',rerrs3,nconds,ndraws);
        prinar2('rerrs4=',rerrs4,nconds,ndraws);


%
%        average errors
%
        rerrs_mean = mean(rerrs,2);
        rerrs2_mean = mean(rerrs2,2);
        rerrs3_mean = mean(rerrs3,2);
        rerrs4_mean = mean(rerrs4,2);
        rerrs5_mean = mean(rerrs5,2);

        prin2('rerrs=',rerrs_mean,nconds);
        prin2('rerrs2=',rerrs2_mean,nconds);
        prin2('rerrs3=',rerrs3_mean,nconds);
        prin2('rerrs4=',rerrs4_mean,nconds);
        prin2('rerrs5=',rerrs5_mean,nconds);


%
%        rank estimates
%
        for i=1:nconds
%
        khat_maxes(i,1) = max(khats(i,:));
        khat_maxes(i,2) = max(khats2(i,:));
        khat_maxes(i,3) = max(khats3(i,:));
        khat_maxes(i,4) = max(khats4(i,:));
        khat_maxes(i,5) = max(khats5(i,:));
%
        khat_mins(i,1) = min(khats(i,:));
        khat_mins(i,2) = min(khats2(i,:));
        khat_mins(i,3) = min(khats3(i,:));
        khat_mins(i,4) = min(khats4(i,:));
        khat_mins(i,5) = min(khats5(i,:));
%
        khat_means(i,1) = mean(khats(i,:));
        khat_means(i,2) = mean(khats2(i,:));
        khat_means(i,3) = mean(khats3(i,:));
        khat_means(i,4) = mean(khats4(i,:));
        khat_means(i,5) = mean(khats5(i,:));
%
        khat_medians(i,1) = median(khats(i,:));
        khat_medians(i,2) = median(khats2(i,:));
        khat_medians(i,3) = median(khats3(i,:));
        khat_medians(i,4) = median(khats4(i,:));
        khat_medians(i,5) = median(khats5(i,:));
    end


        if (0)
%
        prinarf('khats=',khats,nconds,ndraws);
        prinarf('khats2=',khats2,nconds,ndraws);

        prinarf('khat_maxes=',khat_maxes,nconds,4);
        prinarf('khat_mins=',khat_mins,nconds,4);
        prinar2('khat_means=',khat_means,nconds,4);
        prinar2('khat_medians=',khat_medians,nconds,4);

    end


        save('compare_cond.mat','m','n','n0','k','rconds','nconds','rerrs_mean',...
            'rerrs2_mean','rerrs3_mean','rerrs4_mean','rerrs5_mean','khat_means',...
             'khat_mins','khat_maxes','khat_medians','ndraws','time_mc');


        plot_condn();


        end
%
%
%
%
%
        function [rerrs,rerrs2,rerrs3,rerrs4,rerrs5,khats,khats2,khats3,...
            khats4,khats5] = sweep_condnums(rconds,dds,m,n,n0,k,keep,nconds)

        rerrs=zeros(nconds,1);
        rerrs2=zeros(nconds,1);
        rerrs3=zeros(nconds,1);
        rerrs4=zeros(nconds,1);
        rerrs5=zeros(nconds,1);
%
        khats=zeros(nconds,1);
        khats2=zeros(nconds,1);
        khats3=zeros(nconds,1);
        khats4=zeros(nconds,1);
        khats5=zeros(nconds,1);


%%%        prin2('rconds=',rconds,nconds);
%%%        prinf('nconds=',nconds,1);
%%%
%%%        prinf('m=',m,1);
%%%        prinf('n=',n,1);
%%%        prinf('k=',k,1);
%
        gam=m/n;
        gam0=m/n0;

%%%        prin2('dd1=',dds(:,1),10);

%
%        random white noise
%
        epw = randn(m,n);
        epw = epw / sqrt(n);
%
        epw0 = randn(m,n0);
        epw0 = epw0 / sqrt(n0);


%
%        signal singular vectors
%
        [ux,vx] = draw_svecs(m,n,k,dds);


        prin2('epw=',epw,10);
        prin2('epw0=',epw0,10);

        fff = sum(epw.^2,2);
        ggg = sum(epw0.^2,2);

        prin2('fff=',fff,10);
        prin2('ggg=',ggg,10);

%%%        plot(ggg);


%%%        vars0 = make_vars(rconds(1),m);
%%%        vars0 = make_vars(rconds(nconds),m);
%%%        sx = signal_sv(vars0,dds,m,n,n0,k,keep);

%
%        bulk edge
%
        icond=nconds;
        icond=1;
        vars0 = linspace(1,rconds(icond),m)';
        vars0 = vars0 / sum(vars0);

        bvars = 1;
%
        awhts=ones(m,1)/m;
        bwhts=1;
        bedge = mpbdry_edge(vars0,bvars,awhts,bwhts,m,n,gam);
        prin2('bedge=',bedge,1);


        [ell0,cout0,cinn0] = mpbdry_sback(bedge,vars0,bvars,...
            awhts,bwhts,m,n,gam);


%
%        signal singular values
%
        for i=1:k
%
        sx(i) = ell0 + k - i + 1/10;
    end
        prin2('sx=',sx,k);


        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ep,ep0,cov0] = ...
            drawf123(ux,vx,sx,vars0,m,n,k,n0,epw,epw0);

        eee=norm(ep)^2;
        prin2('eee=',eee,1);
        prin2('bedge=',bedge,1);

%%%        prinstop;


        prin2('sx=',sx,k);

        prin2('vars0=',vars0,2);

%%%        prinstop;



        for ijk=1:nconds
%
        rcond=rconds(ijk);

        prin2('rcond=',rcond,1);

        vars2 = make_vars(rcond,m);

        vars = linspace(1,rcond,m)';
        vars = vars / sum(vars);

        prin2('vars=',vars,10);
        prin2('vars2=',vars2,10);

        prinstop;

%%%        sx = signal_sv(vars,dds,m,n,n0,k,keep);

        prin2('vars=',vars,10);
        prin2('vars(m)=',vars(m),1);

        prin2('sx=',sx,k);

        taus = fshr_taus_exact(dds,vars,m,k);
        prin2('taus=',taus,k);

%
%        color the noise and build observed matrix
%
        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ep,ep0,cov0] = ...
            drawf123(ux,vx,sx,vars,m,n,k,n0,epw,epw0);

        prin2('y=',y,10);


%%%        syf_full = svd(yf);
%%%        evals = syf_full.^2;
%%%        fshr_wacompare0(evals,gam,gam0,1);


        [err,err2,err3,err4,err5,xhat,errhat,khat,xhat2,khat2,xhat3,...
            khat3,errhat3,xhat4,khat4,xhat5,khat5] = ...
            allmethods(y,xmat,cov0,ep0,gam0,vars,ep,m,n,k);



        prinf('khat=',khat,1);

        prin2('errhat=',errhat,1);

        prin2('error=',err,1);
        prin2('error, plugin=',err2,1);
        prin2('error, optshrink=',err3,1);
        prin2('error, oracle=',err4,1);

%
%        relative errors
%
        xnorm = norm(xmat,'fro');
        rerr = err / xnorm;
        rerr2 = err2 / xnorm;
        rerr3 = err3 / xnorm;
        rerr4 = err4 / xnorm;
        rerr5 = err5 / xnorm;


        prin2('relative error=',rerr,1);
        prin2('relative error, plugin=',rerr2,1);
        prin2('relative error, optshrink=',rerr3,1);
        prin2('relative error, oracle=',rerr4,1);
        prin2('relative error, plugin2=',rerr5,1);


        rerrs(ijk) = rerr;
        rerrs2(ijk) = rerr2;
        rerrs3(ijk) = rerr3;
        rerrs4(ijk) = rerr4;
        rerrs5(ijk) = rerr5;


        khats(ijk) = khat;
        khats2(ijk) = khat2;
        khats3(ijk) = khat3;
        khats4(ijk) = khat4;
        khats5(ijk) = khat5;


    end

        prin2('rel errs=',rerrs,nconds);
        prin2('rel errs2=',rerrs2,nconds);
        prin2('rel errs3=',rerrs3,nconds);
        prin2('rel errs4=',rerrs4,nconds);
        prin2('rel errs5=',rerrs5,nconds);

%%%        prinstop;
        end
%
%
%
%
%
        function [err,err2,err3,err4,err5,xhat,errhat,khat,xhat2,khat2,xhat3,...
            khat3,errhat3,xhat4,khat4,xhat5,khat5] = ...
            allmethods(y,xmat,cov0,ep0,gam0,vars,ep,m,n,k)

%
%        pseudo-whitening
%
        [xhat,xfhat,errhat,khat,uyf2,vy,tvals,tvals2,taus_hat,sx_hat,couts,...
            cinns,dmu] = fshr_approx(y,cov0,gam0,m,n,k);
        err = norm(xhat - xmat,'fro');
%
%        plug-in covariance
%
        [xhat2,khat2] = plugin44(y,m,n,k,ep0,cov0);
        err2 = norm(xhat2 - xmat,'fro');
%
%        optshrink
%
        [xhat3,khat3,errhat3] = optshrink44(y,m,n,k,vars,ep);
        err3 = norm(xhat3 - xmat,'fro');
%
%        oracle whitening
%
        [xhat4,khat4] = oracle_wshr(y,m,n,k,vars);
        err4 = norm(xhat4 - xmat,'fro');


%
%        plug-in covariance, better rank estimate
%
        prini(0,0);
        [xhat5,khat5] = plugin44(y,m,n,khat,ep0,cov0);
        err5 = norm(xhat5 - xmat,'fro');
        prini(13,0);


        end
%
%
%
%
%
        function [ux,vx] = draw_svecs(m,n,k,dds)
%

        ux = randn(m,k);
        ux=ux.*dds;

%%%        ux = gramschmidt334(ux,m,k);
        ux = fshr_gramschmidt(ux,m,k);
%
        vx = randn(n,k);
%%%        vx = gramschmidt334(vx,n,k);

        vx = fshr_gramschmidt(vx,n,k);

        end
%
%
%
%
%
        function [y,yf,xf,x,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ep,ep0,cov0] = ...
            drawf123(ux,vx,sx,vars,m,n,k,n0,epw,epw0)

%
%        Gaussian noise
%
%%%        ep = diag(sqrt(vars)) * epw;
        ep = fshr_dleft(sqrt(vars),epw,m,n);
%
%        generate signal and spiked matrices
%
        x = ux * diag(sx) * vx';
        y = x + ep;

%
%        independent noise stream
%
%%%        ep0 = diag(sqrt(vars)) * epw0;
        ep0 = fshr_dleft(sqrt(vars),epw0,m,n);

%
%        spiked F-matrix
%
        cov0 = ep0*ep0';
%%%        [cinv0,wmat0,whts] = invdumb222(cov0,m);

        cinv0=inv(cov0);
        wmat0=cinv0^(1/2);
        whts=cov0^(1/2);

        xf = wmat0 * x;
        epf = wmat0 * ep;
        yf = wmat0 * y;

        [uyf,syf,vyf] = fshr_svds(yf,m,n,k);
        [uxf,sxf,vxf] = fshr_svds(xf,m,n,k);


        return;

        evals = svd(ep);
        evals = evals.^2;

        prin2('evals=',evals,10);

        condn = evals(1) / evals(m);
        prin2('condn=',condn,1);

        evals0 = svd(ep0);
        evals0 = evals0.^2;

        prin2('evals0=',evals0,10);

        condn0 = evals0(1) / evals0(m);
        prin2('condn0=',condn0,1);

%%%        prinstop;

        end
%
%
%
%
%
        function sx = signal_sv(vars,dds,m,n,n0,k,keep)

        sx = zeros(k,1);

        taus = fshr_taus_exact(dds,vars,m,k);
        prin2('taus=',taus,k);

        gam = m/n;
        gam0 = m/n0;
        [ycut,xcut] = fshr_cuts(gam,gam0);
        prin2('xcut=',xcut,1);

%%%        prinstop;

        for i=1:keep
%
        sx(i) = xcut*(2 + (keep-i)/keep) / sqrt(taus(i));
        sx(i) = xcut*(1 + (keep-i+1)/keep) / sqrt(taus(i));
%%%        sx(i) = xcut*(2 + k-i) / sqrt(taus(i));
    end

 
        end
%
%
%
%
%
        function vars = make_vars(rcond,m)

        vars=zeros(m,1);

%
%        equispaced variances
%
        for i=1:m
%
%%%        vars(i) = rcond*(m-i)/(m-1) + i/(m-1);
        vars(i) = rcond*(i-1)/(m-1) + (m-i)/(m-1);

    end

%%%        vars = vars / sum(vars);

%
%        rescale variances so that tau is 1
%
        dds = ones(m,1);
        tau = fshr_taus_exact(dds,vars,m,1);

        vars = vars * tau;

        return;

        prin2('tau=',tau,1);
        prin2('vars=',vars,100);

%
%        check that tau is 1
%
        tau = fshr_taus_exact(dds,vars,m,1);
        prin2('tau=',tau,1);

        chk0 = tau-1;
        prin2('chk0=',chk0,1);

        prinstop;
 
        end
%
%
%
%
%
        function [xhat,khat] = oracle_wshr(y,m,n,k,vars)

        y=sqrt(n)*y;
        [xhat,uhat,vhat,t_hat,khat,s_hat,uyw,vyw,syw,ellsw,err_hat,...
            cos_out,cos_inn,cosw_out,cosw_inn,ataus] = wshr_cols(y,m,n,k,vars);
        xhat = xhat / sqrt(n);

        end
%
%
%
%
%
        function [xhat,khat,errhat] = optshrink44(y,m,n,k,vars,ep)

        gam=m/n;

        as=vars;
        bs=1;
        awhts=ones(m,1)/m;
        bwhts=1;
%
        bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam);
%%%        bedge=sqrt(bedge);

        prin2('bedge=',bedge,1);

        ep_norm = norm(ep)^2;
        prin2('ep_norm=',ep_norm,1);


%%%        bedge=0;

        bedge=bedge/n;

        tol=1e-2;
        bedge = bedge*(1+tol);


%%%        prinstop;

        [xhat,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,eee,khat] = ...
           svshr_stiel2(y,m,n,k,bedge);
        errhat=sqrt(sum(eee))*sqrt(n);

%%%        xhat = xhat / sqrt(n);

        prinf('k=',k,1);
        prinf('khat, optshrink=',khat,1);

        end
%
%
%
%
%
        function [xhat,khat] = plugin44(y,m,n,k,ep0,cov0)
%
%        denoise with oracle shrinker and plug-in covariance
%

        [u0,s0,v0] = svd(ep0);
%%%        chk0 = norm(u0*s0*v0' - ep0);

        s0=diag(s0);
        avars = s0.^2;
        avecs = u0;


        gam=m/n;

        if (0)
        prin2('avars=',avars,6);
        prin2('vars_true=',vars_true,6);

        vmean = mean(avars);
        prin2('mean variance=',vmean,1);

        vmin = min(avars);
        prin2('min variance=',vmin,1);
%%%        histogram(avars,10);


        chk0 = norm(avecs * diag(avars) * avecs' - cov0);
        prin2('chk0=',chk0,1);
%%%        prinstop;

        prin2('s0=',s0,20);
        prin2('avars=',avars,20);

        prin2('chk0=',chk0,1);





%%%        avecs = eye(m);


        prin2('vars_true=',vars_true,6);
        prin2('avars=',avars,6);


        prin2('vars_true=',min(vars_true),1);
        prin2('avars=',min(avars),1);

    end

%%%        prinstop;

        y2=y*sqrt(n);
        [xhat,uhat,vhat,t_hat,khat,s_hat,uyw,vyw,syw,ellsw,err_hat,...
            cos_out,cos_inn,cosw_out,cosw_inn,ataus] = wshr_dense(y2,m,n,k,avars,avecs);


        err_hat = sqrt(err_hat);
        xhat=xhat/sqrt(n);


        prinf('khat2=',khat,1);

%%%        prinstop;

        end
%
%
%
%
%
        function compare_wshr()
%
        randn(1,15111);

        rcond = 1000;

        prin2('rcond=',rcond,1);

        m=600;
        n=floor(3*m);
        n0=floor(2*m);

        k=5;
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

%%%        dds=linspace(1,23,m)';
%%%        dds(1:m/2,1) = 3;
%%%        dds(m/2+1:m,1) = 1;
%%%        dds(:,1)=dds(:,1)/norm(dds(:,1)) * sqrt(m);

        prin2('dd1=',dds(:,1),10);

%%%         prinstop;

        taus = fshr_taus_exact(dds,vars,m,k);
        prin2('taus=',taus,k);

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
        sx(i) = xcut*(.01 + (k-i)/k) / sqrt(taus(i));


        sx(i) = (k-i+1)/k + .1;
%%%        sx(i) = xcut*(1.5 + (k-i)) / sqrt(taus(i));
    end
%%%        sx(k)=0;

%%%        sx=zeros(k,1);

        prin2('sx=',sx,k);

%
%        parameter mu (average variance)
%
        dmu=0;
        for i=1:m
%
        dmu = dmu+vars(i)/m;
    end

        prin2('dmu=',dmu,1);

%%%        [y,yf,xf,xmat,epf,whts,wmat0,uyf,syf,vyf,uxf,sxf,vxf,ux,vx,ep,ep0,cov0,epw,epw0] = ...
%%%            fshr_draw(sx,vars,m,n,k,n0,dds);


        [xmat,ux,vx,y,yw,uyw,vyw,syw] = draw_spiked333(sx,vars,m,n,k);

        xmat = xmat / sqrt(n);
        y = y / sqrt(n);

        prin2('vars=',vars,10);
        prin2('sx=',sx,k);

        sall = svd(yw);
%
%%%        whtd_mpcompare(yw/sqrt(n),m,n,1,1);
%%%        plot(sall,'*')

        prin2('xmat=',xmat,10);

%%%        prinstop;

        prin2('xmat=',xmat,10);

%%%        prinstop;

%%%        prin2('y=',y,100);
%%%        prin2('yf=',yf,100);
%%%        prin2('xmat=',xmat,100);
%%%        prin2('whts=',whts,100);
%%%        prin2('epf=',epf,100);
%%%        prin2('epw0=',epw0,100);


%
%        new method
%

        gam0=0;

        [xhat,xfhat,errhat,khat,uyf2,vy,tvals,tvals2,taus_hat,sx_hat,couts,...
            cinns,dmu] = fshr_approx(y,diag(vars),gam0,m,n,k,iorth);


        y2=sqrt(n)*y;
        [xhat2,uhat2,vhat2,t_hat2,khat2,s_hat2,uyw2,vyw2,syw2,ellsw2,err_hat2,...
            cos_out2,cos_inn2,cosw_out2,cosw_inn2,ataus2] = wshr_cols(y2,m,n,k,vars);
        xhat2 = xhat2 / sqrt(n);


        prin2('tvals2=',tvals2,k);
        prin2('t_hat2=',t_hat2,k);


        ddd = tvals2 - t_hat2';
        prin2('ddd=',ddd,k);

%%%        prinstop;

%%%        prin2('vy=',vy,10);
%%%        prin2('vyw2=',vyw2,10);

        prin2('cinns=',cinns,k);
        prin2('cos_inn2=',cos_inn2,k);

        prinf('khat=',khat,1);



        ddd = cinns - cos_inn2';
        prin2('ddd=',ddd,k);


        prin2('taus_hat=',taus_hat,k);
        prin2('ataus2=',ataus2,k);


        ddd = taus_hat-ataus2';
        prin2('ddd=',ddd,k);

%%%        prinstop;

        prin2('xhat=',xhat,100);
        prin2('xfhat=',xfhat,100);
        prin2('tval=',tvals2,k);
        prin2('taus_hat=',taus_hat,k);
        prin2('cws=',cws,k);
        prin2('cinns=',cinns,k);
        prin2('uyf2=',uyf2,100);
        prin2('vy=',vy,100);


        err = norm(xhat - xmat,'fro');

        prin2('errhat=',errhat,1);
        prin2('error=',err,1);

        prin2('xhat2=',xhat2,5);
        prin2('xhat=',xhat,5);


        prin2('errhat=',errhat,1);
        prin2('err_hat2=',err_hat2,1);

        dif = xhat2 - xhat;
        prin2('dif=',dif,10);


        chk0 = norm(dif,'fro') / norm(xhat,'fro');
        prin2('chk0=',chk0,1);


        prin2('errhat=',errhat,1);
        prin2('err_hat2=',sqrt(err_hat2),1);


        chk0 = errhat - sqrt(err_hat2);
        prin2('chk0=',chk0,1);



        prinstop;




        [xhat3,xfhat3,errhat3,khat3,tvals2_3,sx3,cws3,cinns3,uyf2_3] = ...
            fshr_core(yf,whts,m,n,n0,k,taus_hat,dmu);


        ddd = xhat3 - xhat;
        chk0 = norm(ddd,'fro');
        prin2('ddd=',ddd,10);
        prin2('chk0=',chk0,1);

        ddd = xfhat3 - xfhat;
        chk0 = norm(ddd,'fro');
        prin2('chk0=',chk0,1);

        prinf('khat=',khat,1);
        prinf('k=',k,1);

        err = norm(xhat - xmat,'fro');

        prin2('errhat=',errhat,1);
        prin2('error=',err,1);


        prinstop;
        end
%
%
%
%
%
        function [x,ux,vx,y,yw,uyw,vyw,syw] = draw_spiked333(sx,vars,m,n,k)
%
%        signal singular vectors
%
        ux = randn(m,k);
%%%        ux = gramschmi2(ux,m,k);
        ux = fshr_gramschmidt(ux,m,k);
%
        vx = randn(n,k);
%%%        vx = gramschmi2(vx,n,k);


        chk0 = norm(ux'*ux - eye(k),'fro');
        prin2('chk0=',chk0,1);

        chk0 = norm(vx'*vx - eye(k),'fro');
        prin2('chk0=',chk0,1);



%
%        generate Gaussian noise
%
        epw = randn(m,n);
%%%        epw = epw / sqrt(n);
%%%        ep = diag(sqrt(vars)) * epw;

        ep = wshr_dlmult(epw,sqrt(vars),m,n);


%
%        (check variances are correct)
%
        i=2;
        vi = sum(ep(i,:).^2)/n;

        prin2('vi=',vi,1);
        prin2('variance=',vars(i),1);


%%%        prinstop;

%
%        generate signal and spiked matrices
%
        x = ux * diag(sx) * vx';
        y = x + ep;


%


        yw = wshr_dlmult(y,sqrt(1./vars),m,n);
        xw = wshr_dlmult(x,sqrt(1./vars),m,n);

        chk0 = norm(xw+epw - yw,'fro');
        prin2('chk0=',chk0,1);

%
%        SVDs of matrices
%
        [ux,sx,vx] = svshr_svdsmart(x,m,n,k);

        chk0 = norm(ux*diag(sx)*vx' - x,'fro');
        prin2('chk0=',chk0,1);


        [uyw,syw,vyw] = svshr_svdsmart(yw,m,n,k);

        chk0 = norm(yw*vyw - uyw*diag(syw),'fro');
        prin2('chk0=',chk0,1);

        chk0 = norm(yw'*uyw - vyw*diag(syw),'fro');
        prin2('chk0=',chk0,1);



%%%        prinstop;



        end
%
%
%
%
%
        function err = err_oracle(sx,tau,dmu,gam)

        amu=dmu;
        alpha = 1/tau;
        sx2 = sx*sqrt(tau);
%%%        sx2 = sx*tau;

        bmu=1;
        beta=1;
        sig=1;
        [err,rlam,cout,cinn] = whtd_errfmla(sx2,alpha,beta,gam,sig,amu,bmu);

%%%        err=err*sqrt(n);
%%%        prin2('err=',err,1);



        end
%
%
%
%
%
        function tau = tau_approx(ep0,m,n0)
%
        gam0 = m/n0;

        sep0 = svd(ep0);
        sep0 = sep0 / sqrt(n0);
        tau = sum(1./sep0.^2) / m  * (1-gam0);

%%%        tau = 1 / tau / (1-gam0);

        prin2('tau, estimated=',tau,1);


        end
%
%
%
%
%
        function recipr_compare555(evals,gam,ifig)
%
%        plots histogram of evals and overlays plot of density of
%        reciprocal eigenvalues for white noise with variance 1, with aspect
%        ratio gam. ifig is the figure number of plot
%

        [f,x0,x1] = mpreceval(1000,gam);


        nvals=100;
        vals = zeros(1,nvals);
        ts=linspace(x0,x1,nvals);

        for i=1:nvals
%
        [vals(i),bpls,bmin] = mpreceval(ts(i),gam);
    end


        evals = 1./evals;

        figure(ifig);
        subplot(1,2,1);
        hh = histogram(evals,100,'Normalization','pdf');
        hold on; plot(ts,vals,'linewidth',3)

        m=length(evals);

        subplot(1,2,2);
        plot(evals,'*','Color','b')
        hold on;
        plot(0:m+1,x1*ones(m+2,1),'LineWidth',2,'Color','r')
        plot(0:m+1,x0*ones(m+2,1),'LineWidth',2,'Color','r')

        xlim([0,m+1])

        set(figure(ifig),'Position',[500,500,1300,500])

        end
%
%
%
%
%
        function rint = stint_dumb2(x,gam,gam0)
%
%        evaluate Wachter's Stieltjes transform numerically
%
        prin2('x=',x,1);

        [bmin,bmax] = fshr_wachter_lims(gam,gam0);

        fh = @(t) eval_integrand(t,x,gam,gam0);
        rint = integral(fh,bmin,bmax);

        end
%
%
%
%
%
        function fs = eval_integrand(xs,x,gam,gam0)

        n=length(xs);
        fs = zeros(size(xs));
        for i=1:n
%
        fs(i) = fshr_evalwacht(xs(i),gam,gam0);
        fs(i) = fs(i) / (xs(i)-x);
    end

        end
%
%
%
%
%
        function [s,der,sbar,dbar] = wachtst_left(x,gam,gam0)
%
%        evaluate Wachter Stieltjes transform for x on the left
%
        sss1 = 1-gam-(1+gam0)*x;
        sss2 = sqrt(((1-gam0)*x + 1-gam)^2 - 4*x);
        sss3 = 2*x*(gam+x*gam0);

        s = sss1 - sss2;
        s = s / sss3;
%
%        evaluate the derivative
%
        ddd1 = -1-gam0;
        ddd2 = 2*((1-gam0)*x + (1-gam))*(1-gam0) - 4;
        ddd2 = ddd2/sss2/2;
        ddd3 = 2*gam + 4*gam0*x;

        der = sss3*(ddd1-ddd2) - (sss1-sss2)*ddd3;
        der = der / sss3^2;

%
%        evaluate the complementary transform
%
        sbar = gam*s + gam/x - 1/x;
        dbar = gam*der - gam/x^2 + 1/x^2;

        end
%
%
%
%
%
        function [s,der,sbar,dbar] = wachtst_right(x,gam,gam0)
%
%        evaluate Wachter Stieltjes transform for x on the right
%
        sss1 = 1-gam-(1+gam0)*x;
        sss2 = sqrt(((1-gam0)*x + (1-gam))^2 - 4*x);
        sss3 = 2*x*(gam+x*gam0);

        s = sss1 + sss2;
        s = s / sss3;

%
%        evaluate the derivative
%
        ddd1 = -1-gam0;
        ddd2 = 2*((1-gam0)*x + (1-gam))*(1-gam0) - 4;
        ddd2 = ddd2/sss2/2;
        ddd3 = 2*gam + 4*gam0*x;

        der = sss3*(ddd1+ddd2) - (sss1+sss2)*ddd3;
        der = der / sss3^2;

%
%        evaluate the complementary transform
%
        sbar = gam*s + gam/x - 1/x;
        dbar = gam*der - gam/x^2 + 1/x^2;

        end
%
%
%
%
%
        function [f,bpls,bmin] = mpreceval(x,gam)
%
%        reciprocal Marchenko Pastur density, sigma=1
%
        x0 = (1-sqrt(gam))^2;
        x1 = (1+sqrt(gam))^2;

        bpls = 1/x0;
        bmin = 1/x1;

        f = (bpls - x) * (x - bmin);
        f = sqrt(f*x0*x1) / (2*pi*gam) / x^2;

        end
%
%
%
%
%
        function [yval,yder] = spikforw_der(xval,gam,gam0)

        xx=xval^2;
%
%        spike-forward map
%
        ytop = (1 + xx) * (gam + xx);
        ybot = (1-gam0)*xx - gam0;
        yval = ytop / ybot;
%
%        derivative
%
        ytop_der = (1+gam+2*xx);
        ybot_der = 1-gam0;
        yder = (ybot*ytop_der - ytop*ybot_der) / ybot^2;
        yder = yder*2*xval;

%
%        take square root and adjust derivative
%
        yder = yder / sqrt(yval) / 2;
        yval = sqrt(yval);

        end
%
%
%
%
%
        function [res,der] = resolv_left(yeig,gam,gam0)
%
        rr = sqrt(gam+gam0-gam*gam0);
        bb = (1-gam0)*yeig - (1-gam);
        dd = bb^2 - 4*yeig*rr^2;
        res = -2 / (bb-sqrt(dd));

        bder = 1-gam0;
        dder = 2*bb*bder-4*rr^2;
        der = 2*(bder - dder/sqrt(dd)/2) / (bb-sqrt(dd))^2;

        return;
%
%        alternate formula
%

        res2 = res;
        res = ((1-gam0)*yeig + 1 - gam)^2 - 4*yeig;
        res = 1-gam-yeig*(1-gam0) - sqrt(res);
        res = res / 2 / (gam + gam0 - gam*gam0) / yeig;

        chk0 = res-res2;
        prin2('res=',res,1);
        prin2('res2=',res2,1);
        prin2('chk0=',chk0,1);


        prinstop;

        end
%
%
%
%
%
        function [res,der] = resolv_right(yeig,gam,gam0)
%
        rr = sqrt(gam+gam0-gam*gam0);
        bb = (1-gam0)*yeig - (1-gam);
        dd = bb^2 - 4*yeig*rr^2;
        res = -2 / (bb+sqrt(dd));

        bder = 1-gam0;
        dder = 2*bb*bder-4*rr^2;
        der = 2*(bder + dder/sqrt(dd)/2) / (bb+sqrt(dd))^2;


        return;

%
%        alternate formula
%

        res2 = res;
        res = ((1-gam0)*yeig + 1 - gam)^2 - 4*yeig;
        res = 1-gam-yeig*(1-gam0) + sqrt(res);
        res = res / 2 / (gam + gam0 - gam*gam0) / yeig;

        chk0 = res-res2;
        prin2('res=',res,1);
        prin2('res2=',res2,1);
        prin2('chk0, inside resolv=',chk0,1);

        prinstop;
%
        end
