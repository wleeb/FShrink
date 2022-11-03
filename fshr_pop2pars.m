        function [tval,eta,err,syf,cinn,cout,cc,cw,upone,dupone,uptwo] = ...
            fshr_pop2pars(sx,dmu,gam,gam0,tau)

        cout=0;
        cinn=0;
        cc=0;
        cw=0;
%
        tval=0;
        eta=0;
%
        err=0;
        unorm=0;
        upone=0;
        dupone=0;
        uptwo=0;

        syf=0;

        if (tau <= 0)
%
        return;
    end

        xcut_tau = fshr_xcut_tau(gam,gam0,tau);
        if (sx <= xcut_tau)
%
        return;
    end

%
%        empirical singular value and eigenvalue
%
        syf = fshr_spikforw(sx,gam,gam0,tau);
        yy=syf^2;

%
%        key transforms
%
        [upone,dupone] = fshr_upone_fmla(yy,gam,gam0);
        uptwo = fshr_uptwo_fmla(yy,gam,gam0);

        [vpsi,dpsi] = fshr_evalpsi(yy,gam,gam0);
        [sbar,sbder] = fshr_wachtbarst(yy,gam,gam0);

        [ss,sder] = fshr_wachtstiel_right(yy,gam,gam0);
        [res,der] = fshr_resolv(yy,gam,gam0);

%
%        vector norm
%
        t1 = -gam*abs(res)*(upone + yy*dupone);
        t2 = yy*abs(sbar)*(uptwo-ss^2);
        efun = t1 + t2;
        vphi = yy*abs(sbar)*ss^2;
        unorm=sqrt((vphi/tau + efun*dmu) / abs(dpsi));

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
%        optimal singular value and square error
%
        err = sx^2*(1 - cinn^2*cw^2 / unorm^2);
        eta = sx*cinn*cw / unorm^2;
        tval = sx*cinn*cout;
%%%        tval = eta*unorm;


        end
