        function [tval,eta,err,sx,cinn,cout,cc,cw,upone,dupone,uptwo,tau] = ...
            fshr_emp2pars(syf,unorm,dmu,gam,gam0)

        cout=0;
        cinn=0;
        cc=0;
        cw=0;
%
        sx=0;
        tval=0;
        eta=0;
%
        err=0;
        upone=0;
        dupone=0;
        uptwo=0;
        tau=0;
%
        [ycut,xcut] = fshr_cuts(gam,gam0);
        if (syf <= ycut)
%
        return;
    end

        yy=syf^2;

%
%        key transforms
%
        [vpsi,dpsi] = fshr_evalpsi(yy,gam,gam0);
        [sbar,sbder] = fshr_wachtbarst(yy,gam,gam0);

        [ss,sder] = fshr_wachtstiel_right(yy,gam,gam0);
        [res,der] = fshr_resolv(yy,gam,gam0);

        [upone,dupone] = fshr_upone_fmla(yy,gam,gam0);
        uptwo = fshr_uptwo_fmla(yy,gam,gam0);
%

%
%        formula for tau
%
        t1 = -gam*abs(res)*(upone + yy*dupone);
        t2 = yy*abs(sbar)*(uptwo-ss^2);
        efun = t1 + t2;

%%%        prin2('t2=',t2,1);
%%%        prin2('t1=',t1,1);
%%%        prin2('efun=',efun,1);

        vphi = yy*abs(sbar)*ss^2;

        tau = vphi / (abs(dpsi)*unorm^2 -  efun*dmu);
%
%        exit if tau is negative
%
        if (tau <= 0)
%
        return;
    end

%
%        inner products
%
        cinn = sbar*vpsi / dpsi;
        cinn = sqrt(cinn);
%
        cc = tau*res*vpsi/dpsi;
        cc=sqrt(cc);
%
        cw = yy * sbar * ss^2 / tau / dpsi;
        cw = sqrt(cw);
        cout = cw / unorm;

%
%        singular values and squared error
%
        sx = fshr_spikeback(syf,gam,gam0,tau);
        err = sx^2*(1 - cinn^2*cw^2 / unorm^2);
        eta = sx*cinn*cw / unorm^2;
        tval = sx*cinn*cout;
%%%        tval = eta*unorm;


        end
