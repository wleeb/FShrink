        function main
        prini(13,1);
%

%%%        test444();
%%%        prinstop;

        randn(1,1000);

        m=3000;
        n=2000;
%%%        n=floor(m/100);

        prinf('m=',m,1);
        prinf('n=',n,1);

        gam = m / n;

        as = rand(1,m);
        bs = rand(1,n);

        awhts=ones(1,m);
        awhts = awhts / sum(awhts);
        bwhts=ones(1,n);
        bwhts = bwhts / sum(bwhts);

        rr=1000;
        for i=1:m
%
%%%        as(i) = rr*i/m;
    end

        prin2('as=',as,m);

        difa = sum(awhts) - 1;
        difb = sum(bwhts) - 1;
        prin2('difa=',difa,1);
        prin2('difb=',difb,1);

%%%        prinstop;

   
        bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam);

        prin2('bedge=',bedge,1);


        krank=3;


        ells = bedge + rand(1,krank);
        ells = sort(ells,'descend');
        [y,x,ep,u,v] = genspiket(ells,m,n,krank,as,awhts,bs,bwhts);


        [uy,sy,vy] = svdsmart(y,m,n,krank);


        chk0 = norm(y*vy - uy*diag(sy),'fro') / norm(sy,'fro');

        prin2('chk0=',chk0,1);



%
%        compare asyptotic and observed edge
%
        bedge_obs = norm(ep)^2;
        err = (bedge_obs - bedge) / bedge;

        prin2('edge,observed=',bedge_obs,1);
        prin2('edge,exact=',bedge,1);
        prin2('relative error=',err,1);



%
%        compare spiked model parameters
%
        i=2;
        [rlam,cout,cinn] = mpbdry_sforw(ells(i),as,bs,...
            awhts,bwhts,m,n,gam);
        rlam_obs = sy(i)^2;
        err = (rlam_obs - rlam) / rlam;

        prin2('lambda,observed=',rlam_obs,1);
        prin2('lambda,exact=',rlam,1);
        prin2('relative error=',err,1);


        prinstop;


        cout_obs = abs(sum(u(:,i).*uy(:,i)));
        err = (cout_obs - cout) / cout;


        prin2('cosine,observed=',cout_obs,1);
        prin2('cosine,exact=',cout,1);
        prin2('relative error=',err,1);


        cinn_obs = abs(sum(v(:,i).*vy(:,i)));
        err = (cout_obs - cout) / cout;


        prin2('cosine,observed=',cinn_obs,1);
        prin2('cosine,exact=',cinn,1);
        prin2('relative error=',err,1);

%%%        prinstop;

%
%        check spike backwards mapping
%
        [ell2,cout2,cinn2] = mpbdry_sback(rlam,as,bs,...
            awhts,bwhts,m,n,gam);

        chk0 = abs(ell2 - ells(i));
        chk0 = chk0+abs(cout2 - cout);
        chk0 = chk0+abs(cinn2 - cinn);


        prin2('chk0, backwards mapping=',chk0,1);




        prinstop;
        end
%
%
%
%
%
        function test444()

        m=100;
        n=2000;

        prinf('m=',m,1);
        prinf('n=',n,1);

        gam = m / n;

        as = rand(1,m);
        bs = rand(1,n);

        awhts=ones(1,m);
        awhts = awhts / sum(awhts);
        bwhts=ones(1,n);
        bwhts = bwhts / sum(bwhts);


        difa = sum(awhts) - 1;
        difb = sum(bwhts) - 1;
        prin2('difa=',difa,1);
        prin2('difb=',difb,1);

%%%        prinstop;

%%%        prin2('as=',as,10);
%%%        prin2('bs=',bs,10);


%%%        as=ones(1,m);
%%%        bs=ones(1,n);

        [ell,bedge,err] = mpbdry_thresh(as,bs,awhts,bwhts,m,n,gam);
        ell2 = mpbdry_sback(bedge,as,bs,awhts,bwhts,m,n,gam);

        dif = (ell-ell2)/ell2;
        prin2('dif=',dif,1);

        prin2('fval, should be zero=',err,1);

        prinstop;
        end
%
%
%
%
%
        function [s,sder,evar,eder,gval,gder,gderl] = stielauto(rlam,...
            as,bs,awhts,bwhts,m,n,gam)
%
%         evaluates the Stieltjes transform using fzero (MATLAB root-finding
%         method).
%

        if (size(as,1) ~= 1)
%
        as = as';
    end
%
        if (size(bs,1) ~= 1)
%
        bs = bs';
    end
%
        if (size(awhts,1) ~= 1)
%
        awhts = awhts';
    end
%
        if (size(bwhts,1) ~= 1)
%
        bwhts = bwhts';
    end

        deps = mpbdry_machzero();
%%%        tol=100*deps;
        tol=sqrt(deps)/10;

        feval = @(x) mpbdry_evalf(x,rlam,as,bs,...
            awhts,bwhts,m,n,gam,0);
        evar = fzero(feval,[-1/gam/max(bs)+tol,-tol]);

        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(evar,rlam,as,bs,awhts,...
            bwhts,m,n,gam);

        eder = -fdlam / fder;
        gderl = gder * eder;
%
%        evaluate the stieltjes transform and derivative
%
        s = sum(1 ./ (gval*as-rlam).*awhts);
        sder = -sum(awhts.*(gderl*as - 1) ./ (gval*as - rlam).^2);


        end
%
%
%
%
%
        function [y,x,ep,u,v] = genspiket(ells,m,n,k,as,awhts,bs,bwhts)
%
        ep = randn(m,n) / sqrt(n);
        ep = diag(sqrt(as)) * ep * diag(sqrt(bs));

        u = randn(m,k);
        v = randn(n,k);


        u = gramschmidt(u,m,k);
        v = gramschmidt(v,n,k);

        chk0 = norm(u'*u - eye(k),'fro');
        chk0 = chk0 + norm(v'*v - eye(k),'fro');

        prin2('chk0=',chk0,1);
%%%        prinstop;


        x = u * diag(sqrt(ells)) * v';

        xyz = svds(x,k)';
        prin2('xyz=',xyz,k);
        prin2('xyz=',sqrt(ells),k);
%%%        prinstop;

        y = x + ep; 

        end
%
%
%
%
%
        function u = gramschmidt(uin,m,k)
%
        u = uin;

        u(:,1) = u(:,1) / norm(u(:,1));

%%%        dif = norm(u(:,1))-1;
%%%        prin2('dif=',dif,1);
%%%        prinstop;

        for i=2:k
%
        for ijk = 1:2
%
        for j=1:i-1
%
        pij = sum(u(:,i).*u(:,j));
        u(:,i) = u(:,i) - pij*u(:,j);
    end
        u(:,i) = u(:,i) / norm(u(:,i));
    end

    end

        chk0 = norm(u'*u - eye(k),'fro');
        prin2('chk0=',chk0,1);
%%%        prinstop;

        end
%
%
%
%
%
        function stra = sbar2stra(sbar,z,gam)
%
%        given sbar(z), computes stra(z) (Stieltjes transform)
%
        stra = (sbar + (1-gam)/z) / gam;

        end
%
%
%
%
%
        function stra = stieltjes_mp(z,gam,sig)
%
%        evaluate the Stieltjes transform of standard MP law
%
        b = sig^2*(1-gam) - z;
        disc = (z - sig^2 - gam*sig^2)^2 - 4*gam*sig^4;
        a = 2*gam*z*sig^2;

        stra = (b + sqrt(disc)) / a;

        end
%
%
%
%
%
        function [stra,stra_der] = svals2integrs(s,m,n,rlam);
%
%        computes empirical estimates of the stieljes and D-transforms 
%        at the value rlam, with empirical singular values s (already
%        normalized by dimension).
%
        stra = 0;
        stra_der = 0;

        lmin = min(m,n);
%
        for i=1:lmin
%
        stra = stra + 1 / (s(i)^2 - rlam);
        stra_der = stra_der + 1 / (rlam - s(i)^2)^2;        
    end

        stra = stra/lmin;
        stra_der = stra_der/lmin;

        if (n <= m)
%
        stra = (n/m)*stra - (m-n)/m/rlam;
        stra_der = (n/m)*stra_der + (m-n)/m/rlam^2;
    end

        end
%
%
%
%
%
        function [u,s,v] = svdsmart(a,m,n,k)
%
        if (m/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
    end

        if (m/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
    end

        end
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        This code contains four user-callable functions for analyzing random
%        matrices of the form
%
%                      N = A^{-1/2} G B^{-1/2} / sqrt(n)              (1)
%
%        where A and B are two user-specified sequences of positive weights
%        of lengths m and n, respectively, each with an associated probability
%        distribution. Here, G is an m-by-n matrix of iid N(0,1) entries.
%
%     mpbdry_stiel - evaluates the Stieltjes transform of the LSD of N at a
%        specified value lambda
%     mpbdry_edge - evaluates the asymptotic operator norm squared of the random 
%        matrix N
%     mpbdry_thresh - evaluates the asymptotic operator norm squared of the random 
%        matrix N and the minimum detectable signal eigenvalue
%     mpbdry_sback - evaluates population parameters for spiked matrix 
%        model; from observed value to population value
%     mpbdry_sforw - evaluates population parameters for spiked matrix 
%        model; from population value to observed value
%
%        Additional routines evaluate higher-order derivatives of certain
%        intermediate functions, which are useful for theoretical purposes.
%        Among these routines are:
%
%     mpbdry_evalqder2 - evaluate first and second derivatives of Q
%     mpbdry_fmin - evaluates the minimum of F
%     mpbdry_fminder - evaluates the derivative of the minimum of F
%     mpbdry_evalfmix - evaluates second and mixed derivatives of F
%     mpbdry_evalfg3 - evaluates all derivatives up to third order of both
%        F and G
%     mpbdry_evalfg3_slow - same as mpbdry_evalfg3, but not vectorized
%
%        The methods for these codes are described in the paper
%        ``Rapid evaluation of the spectral signal detection threshold and
%            Stieltjes transform'',
%        by W. Leeb.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
        function [ell,bedge,err] = mpbdry_thresh(as,bs,awhts,bwhts,m,n,gam)
%
%                            description:
%
%   This code evaluates model parameters for random matrices of the form
%
%                      N = A^{-1/2} G B^{-1/2} / sqrt(n)              (1)
%
%   where G has iid Gaussian entries of variance 1. Returns top eigenvalue
%   of YY^T and the minimum detectable signal eigenvalue. The code is
%   based on the observation that the function F(lambda^*,e) (as function of e)
%   appears to have a double root when lambda^* is the top eigenvalue of YY^T;
%   however, this fact has not been proven in general.
%
%
%                           input parameters:
%
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   ell - estimated eigenvalue of XX^T
%   bedge - estimated eigenvalue of YY^T
%   err - the value of F(bedge,evar); this should be zero
%
%
%
%
%        . . . check input dimensions
%
        if (size(as,1) ~= 1)
%
        as = as';
    end
%
        if (size(bs,1) ~= 1)
%
        bs = bs';
    end
%
        if (size(awhts,1) ~= 1)
%
        awhts = awhts';
    end
%
        if (size(bwhts,1) ~= 1)
%
        bwhts = bwhts';
    end

%
%        . . . evaluate the edge of empirical covariance
%
        bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam);

%
%        evaluate the minimum of F(bedge,e) (as function of e)
%
        evar = mpbdry_fmin(bedge,as,bs,awhts,bwhts,m,n,gam);
        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(evar,bedge,as,bs,awhts,...
            bwhts,m,n,gam);
        eder = -fdlam / fder;
        gderl = gder * eder;
%
%        evaluate the stieltjes transform and derivative
%
        s = sum(1 ./ (gval*as-bedge).*awhts);
        sder = -sum(awhts.*(gderl*as - 1) ./ (gval*as - bedge).^2);
        [sbar,sbar_der] = mpbdry_stra2sbar(s,sder,bedge,gam);
        d = s * sbar * bedge;
%%%        dder = sder*sbar*bedge + s*sbar_der*bedge + s*sbar;

        ell = 1 / d;
%
%        evaluate F(bedge,evar) (should be zero)
%
        ifder=0;
        [err,fder] = mpbdry_evalf(evar,bedge,as,bs,awhts,bwhts,m,n,gam,ifder);


        end
%
%
%
%
%
        function [ell,cout,cinn] = mpbdry_sback(rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%                            description:
%
%   This code evaluates model parameters for random matrices of the form
%
%                          Y = X + N                               (1)
%
%   where X is a low rank matrix and N = A^{1/2} G B^{1/2} has separable
%   variance profile. Returns eigenvalue of XX^T and angles between
%   singular vectors of X and Y.
%
%                           input parameters:
%
%   rlam - spiked eigenvalue of YY^T
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   ell - estimated eigenvalue of XX^T
%   cout - estimated cosine between left singular vectors of X and Y
%   cinn - estimated cosine between right singular vectors of X and Y
%
%
%        . . . Stieltjes transforms and D transform
%        (derivative of D wrt rlam, NOT sqrt(rlam))
%
        [s,sder] = mpbdry_stiel(rlam,as,bs,awhts,bwhts,m,n,gam);
        [sbar,sbar_der] = mpbdry_stra2sbar(s,sder,rlam,gam);
        d = s * sbar * rlam;
        dder = sder*sbar*rlam + s*sbar_der*rlam + s*sbar;
%
%        spiked model parameters
%
        ell = 1 / d;
        cout = sqrt(s / (dder * ell));
        cinn = sqrt(sbar / (dder * ell));

        end
%
%
%
%
%
        function [rlam,cout,cinn] = mpbdry_sforw(ell,as,bs,awhts,bwhts,m,n,gam)
%
%                            description:
%
%   This code evaluates model parameters for random matrices of the form
%
%                          Y = X + N                               (1)
%
%   where X is a low rank matrix and N = A^{1/2} G B^{1/2} has separable
%   variance profile. Returns asymptotic eigenvalue of YY^T and angles
%   between singular vectors of X and Y.
%
%                           input parameters:
%
%   ell - eigenvalue of XX^T
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   rlam - estimated eigenvalue of YY^T
%   cout - estimated cosine between left singular vectors of X and Y
%   cinn - estimated cosine between right singular vectors of X and Y
%
%
        mmm=1000;

        rlams = zeros(1,mmm);
        vals = zeros(1,mmm);
        ders = zeros(1,mmm);
%
        bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam);
%
%        solve by Newton
%
        zero_mach = mpbdry_machzero();
        tol=sqrt(zero_mach)/10;
        kstop=0;

        rlams(1)=bedge;
        nsteps=mmm;
%
        for i=1:mmm
%
        [val,vder] = mpbdry_dreciproc(rlams(i),as,bs,awhts,bwhts,m,n,gam);
        vals(i)=val - ell;
        ders(i)=vder;
%
        rlams(i+1) = rlams(i) - vals(i) / ders(i);
%
%        check convergence
%
        rr = (rlams(i+1)-rlams(i)) / rlams(i);
        rr=abs(rr);
%%%        prin2('rr=',rr,1);
        if (rr < tol)
%
        kstop=kstop+1;
    end
        if (kstop==2)
%
        nsteps=i;
        break;
    end
    end

        rlam = rlams(nsteps);
        [ell2,cout,cinn] = mpbdry_sback(rlam,as,bs,awhts,bwhts,m,n,gam);
%%%        chk0 = ell2 - ell

        end
%
%
%
%
%
        function [val,vder] = mpbdry_dreciproc(rlam,as,bs,awhts,bwhts,m,n,gam)
%
        [s,sder] = mpbdry_stiel(rlam,as,bs,awhts,bwhts,m,n,gam);
        [sbar,sbar_der] = mpbdry_stra2sbar(s,sder,rlam,gam);
        d = s * sbar * rlam;
%
%        . . . derivative of D wrt rlam, NOT sqrt(rlam)
%
        dder = sder*sbar*rlam + s*sbar_der*rlam + s*sbar;
%
%        1/D, and its derivative
%
        val = 1/d;
        vder = -dder/d^2;

        end
%
%
%
%
%
        function [sbar,sbar_der] = mpbdry_stra2sbar(stra,stra_der,z,gam)
%
%        given s(z) and s'(z), computes sbar(z) and sbar'(z)
%
        sbar = gam*stra - (1-gam)/z;
        sbar_der = gam*stra_der + (1-gam) / z^2;

        end
%
%
%
%
%
        function [s,sder] = mpbdry_stiel(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%
%                            description:
%
%   This code evaluates the Stieltjes transform and derivative of the
%   limiting spectral distribution for a random matrix of the form 
%   N = A^{1/2} G B^{1/2}, where G is k-by-l with iid entries and the 
%   spectrum of A and B are drawn from specified discrete distributions.
%
%   This function's purpose is just to check that the input vectors are
%   correctly oriented (necessary for vectorization) and to call
%   mpbdry_stiel0, which does all the actual work.
%
%                           input parameters:
%
%   rlam - the value at which we evaluate s and its derivative
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   s - the Stieltjes transform at rlam
%   sder - the derivative of the Stieltjes transform at rlam
%   evar - the value of e(lambda)
%   eder - the value of e'(lambda)
%   gval - the value of G(lambda)
%   gder - derivative of G wrt e
%   gderl - derivative of G wrt lambda 
%
%

        if (size(as,1) ~= 1)
%
        as = as';
    end
%
        if (size(bs,1) ~= 1)
%
        bs = bs';
    end
%
        if (size(awhts,1) ~= 1)
%
        awhts = awhts';
    end
%
        if (size(bwhts,1) ~= 1)
%
        bwhts = bwhts';
    end

        [s,sder,evar,eder,gval,gder,gderl] = mpbdry_stiel0(rlam,...
            as,bs,awhts,bwhts,m,n,gam);

        end
%
%
%
%
%
        function [s,sder,evar,eder,gval,gder,gderl] = mpbdry_stiel0(rlam,...
            as,bs,awhts,bwhts,m,n,gam)

        x0=0;
        [xs,vals,nsteps] = mpbdry_rootf(x0,rlam,as,bs,awhts,bwhts,m,n,gam);
        evar = xs(nsteps);

        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(evar,rlam,as,bs,awhts,...
            bwhts,m,n,gam);

        eder = -fdlam / fder;
        gderl = gder * eder;

%
%        evaluate the stieltjes transform and derivative
%
        s = sum(1 ./ (gval*as-rlam).*awhts);
        sder = -sum(awhts.*(gderl*as - 1) ./ (gval*as - rlam).^2);


        return

%
%        evaluate them the slow way
%
        s = 0;
        for i=1:m
%
        s = s + 1 / (as(i)*gval - rlam) * awhts(i);
    end

        sder = 0;
        for i=1:m
%
        sder = sder - (as(i)*gderl - 1) / (as(i)*gval - rlam)^2 * awhts(i);
    end

        end
%
%
%
%
%
        function [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(evar,...
            rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates F(lambda,e), its first partial derivatives with
%        respect to e; its derivative wrt lambda; and G(e) and its 
%        first derivatives
%
%
%        . . . G and its derivative (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);

%
%        derivative of F wrt lambda
%
        fdlam = -sum(as .* awhts./ (as.*gval - rlam).^2);

%
%        F and its first derivative (wrt e)
%
        fval = evar - sum(as.*awhts ./ (as*gval - rlam));
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

        end
%
%
%
%
%
        function [xs,fvals,nsteps] = mpbdry_rootf(x0,rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%        Finds the root of F_lambda(e), as a function of e, by Newton
%
        nsteps = 1000;
        xs = zeros(1,nsteps);
        fvals = zeros(1,nsteps);
        fders = zeros(1,nsteps);
%
%        initialize at x0
%
        xs(1)=x0;

        deps = mpbdry_machzero();
%%%        tol=100*deps;
        tol=sqrt(deps)/10;

        kstop=0;
        for i=1:nsteps
%
        ifder=1;
        [fvals(i),fders(i)] = mpbdry_evalf(xs(i),rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder);
        xs(i+1) = xs(i) - fvals(i) / fders(i);

        if (abs(fvals(i)) < tol)
%
        kstop = kstop+1;
    end
        if (kstop == 2)
%
        nsteps = i;
        break;
    end
    end
%
        fvals = fvals(1:nsteps);
        xs = xs(1:nsteps);

        end
%
%
%
%
%
        function [fval,fder] = mpbdry_evalf(evar,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder)
%
%        Evaluates F(lambda,e), its first partial derivative with
%        respect to e
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        fval = evar - sum(as.*awhts ./ (as*gval - rlam));

        if (ifder == 0)
%
        fder=-10000;
        return;
    end

        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

        end
%
%
%
%
%
        function bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam)
%
%
%                            description:
%
%   This code computes the right boundary of the limiting spectral
%   distribution for a random matrix of the form N = A^{1/2} G B^{1/2},
%   where G is k-by-l with iid entries and the spectrum of A and B
%   are drawn from specified discrete distributions.
%
%   This function's purpose is just to check that the input vectors are
%   correctly oriented (necessary for vectorization) and to call
%   mpbdry_edge0, which does all the actual work.
%
%                           input parameters:
%
%   as - 1-by-m vector with spectrum of A
%   bs - 1-by-n vector with spectrum of B
%   awhts - 1-by-m vector with probabilities over as
%   bwhts - 1-by-n vector with probabilities over bs
%   m,n - the lengths of as and bs, respectively
%   gam - the aspect ratio of N (not necessarily m/n)
%
%                          output parameters:
%
%   bedge - the boundary of the limiting spectral distribution
%
%
        if (size(as,1) ~= 1)
%
        as = as';
    end
%
        if (size(bs,1) ~= 1)
%
        bs = bs';
    end
%
        if (size(awhts,1) ~= 1)
%
        awhts = awhts';
    end
%
        if (size(bwhts,1) ~= 1)
%
        bwhts = bwhts';
    end

        bedge = mpbdry_edge0(as,bs,awhts,bwhts,m,n,gam);

        end
%
%
%
%
%
        function bedge = mpbdry_edge0(as,bs,awhts,bwhts,m,n,gam)
%
        nsteps=1000;
        rlams = zeros(1,nsteps);
        qvals = zeros(1,nsteps);
        qders = zeros(1,nsteps);

%
%        initialize rlam to left of root (where Q>0)
%
        rlam0 = max(as)*max(bs) * (1+sqrt(gam))^2;
%%%        rlam0 = max(as)*max(bs) + 1;
        for i=1:1000
%
        qval = mpbdry_evalq(rlam0,as,bs,awhts,bwhts,m,n,gam);
        if (qval > 0)
%
        break;
    end
        rlam0 = rlam0/2;
    end

%
%        use Newton to find the root of Q
%
        rlams(1)=rlam0;

        kstop = 0;
        deps = mpbdry_machzero();
        tol = sqrt(deps)/100;

        for i=1:nsteps
%
        [qvals(i),qders(i)] = mpbdry_evalqder(rlams(i),as,bs,awhts,bwhts,m,n,gam);
        rlams(i+1) = rlams(i) - qvals(i) / qders(i);
%
        if (abs(qvals(i)) < tol)
%
        kstop=kstop+1;
    end
        if (kstop==2)
%
        nsteps=i;
        break;
    end
    end

        rlams=rlams(1:nsteps);
        qvals=qvals(1:nsteps);
        qders=qders(1:nsteps);

        bedge = rlams(nsteps);

        end
%
%
%
%
%
        function qval = mpbdry_evalq(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates Q(lambda) = F(lambda,t(lambda)), where t(lambda) minimizes
%        F(lambda,e) (as a function of e); and its first derivative with
%        respect to lambda
%
        val = mpbdry_fmin(rlam,as,bs,awhts,bwhts,m,n,gam);

        ifder=0;
        [fval,fder] = mpbdry_evalf(val,rlam,as,bs,awhts,bwhts,m,n,gam,ifder);
        qval = fval;

        end
%
%
%
%
%
        function [qval,qder] = mpbdry_evalqder(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates Q(lambda) = F(lambda,t(lambda)), where t(lambda) minimizes
%        F(lambda,e) (as a function of e); and its first derivative with
%        respect to lambda
%
        val = mpbdry_fmin(rlam,as,bs,awhts,bwhts,m,n,gam);
        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(val,rlam,as,bs,...
            awhts,bwhts,m,n,gam);

        qval = fval;
        qder = fdlam;

        end
%
%
%
%
%
        function val = mpbdry_fmin(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates the minimum of value of F(lambda,e) (for fixed rlam)
%
        nsteps = 1000;
        evars = zeros(1,nsteps);
        fders = zeros(1,nsteps);
        fders2 = zeros(1,nsteps);

%
%        initialize to the left of the minimizer (where fder < 0)
%
        emin = mpbdry_leftend(rlam,as,bs,awhts,bwhts,m,n,gam);
        einit=emin+abs(emin);

        ifder2=0;
        for i=1:1000
%
        [fder,fder2] = mpbdry_evalfder2(einit,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder2);

        if (fder <= 0)
%
        break;
    end
        einit = (einit + emin)/2;
    end

%
%        find the minimum using Newton
%
        evars(1) = einit;
        ifder2=1;
        kstop=0;
        deps = mpbdry_machzero();
        tol=sqrt(deps)/100;
%
        for i=1:nsteps
%
        [fders(i),fders2(i)] = mpbdry_evalfder2(evars(i),rlam,as,bs,awhts,...
            bwhts,m,n,gam,ifder2);
        evars(i+1) = evars(i) - fders(i) / fders2(i);
%
        if (abs(fders(i)) < tol)
%
        kstop = kstop+1;
    end
        if (kstop == 2)
%
        nsteps = i;
        break;
    end
    end

        fders=fders(1:nsteps);
        fders2=fders2(1:nsteps);
        evars=evars(1:nsteps);

        val = evars(nsteps);


        end
%
%
%
%
%
        function [fder,fder2] = mpbdry_evalfder2(evar,rlam,as,bs,awhts,...
            bwhts,m,n,gam,ifder2)
%
%        Evaluates the first and second partial derivatives of F(lambda,e) with
%        respect to e; second derivative is returned only if ifder2 is not 0
%
%
%        . . . G and its derivative (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);

%
%        first derivative of F (wrt e)
%
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

        if (ifder2 == 0)
%
        fder2=-1000;
        return
    end
%
%        second derivative of F (wrt e)
%
        term1 = sum(bs.^3 .* bwhts ./ (1+gam*evar*bs).^3);
        term2 = sum(awhts .* as.^2 ./ (as*gval - rlam).^2);
        term3 = sum(awhts .* as.^3 ./ (as*gval - rlam).^3);

        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;

        end
%
%
%
%
%
        function eleft = mpbdry_leftend(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%
%        Computes the left endpoint of the interval I_lambda; that is,
%        the root of the function G(e) - max(as) / lambda on the interval 
%        J with e > -1/(gamma * max(bs)).
%
        bmax = max(bs);
        const = min(rlam ./ as);
%
        nsteps=1000;
        vals = zeros(1,nsteps);
        es = zeros(1,nsteps);

%
%        find starting location to left of the root by bisection
%
        emin = -1 / (gam*bmax);
        einit = emin + abs(emin);

        for i=1:1000
%
        ifder=0;
        [gval,gder] = mpbdry_evalg(einit,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder);
        val = gval - const;
%
        if (val >= 0)
%
        break;
    end

        einit = (einit + emin) / 2;
    end

%
%        find root by Newton's method
%
        es(1) = einit;
%
        deps=mpbdry_machzero();
        dzero=sqrt(deps)/100;
        kbreak=0;
%
        for ijk=1:nsteps
%
        ifder=1;
        [gval,gder] = mpbdry_evalg(es(ijk),rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder);
        val = gval - const;
        der = gder;

        vals(ijk) = val;
        es(ijk+1) = es(ijk) - val / der;
%
%        check if converged
%
        if (abs(val)/const <= dzero)
%
        kbreak=kbreak+1;
    end
        if (kbreak == 2)
%
        nsteps=ijk;
        break;
    end

    end

        vals = vals(1:nsteps) / const;
        eleft = es(nsteps);


        end
%
%
%
%
%
        function [gval,gder] = mpbdry_evalg(evar,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder)
%
%        Evaluates G(e) and its first derivative
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));

        if(ifder==0)
%
        gder=-10000;
        return;
    end

        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);

        end
%
%
%
%
%
        function zero_mach = mpbdry_machzero()
%
        zero_mach=1000;

        d = 2.1d0;
        dd=2;
%%%        d=single(d);
%%%        dd=single(dd);

        for i=1:100000
%
        dd=dd / 2;
        dp = dd+d;
        dif = dp - d;

        if (dif == 0)
%
        zero_mach = dd;
        break;
    end
    end

        end
%
%
%
%
%
        function [qval,qder,qder2] = mpbdry_evalqder2(rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%        Evaluates Q(lambda) = F(lambda,t(lambda)), where t(lambda) minimizes
%        F(lambda,e) (as a function of e); and its first two derivatives with
%        respect to lambda
%
        [val,der] = mpbdry_fminder(rlam,as,bs,awhts,bwhts,m,n,gam);
        [fval,fder,gval,gder,fdlam] = mpbdry_evalfg(val,rlam,as,bs,...
            awhts,bwhts,m,n,gam);

        qval = fval;
        qder = fdlam;
%
        qder2 = 2*(as.*awhts.*(as*gder*der-1) ./ (gval*as - rlam).^3);
        qder2=sum(qder2);

        return
%
%        compute qder2 the slow way:
%
        qder2=0;
        for i=1:m
%
        qder2 = qder2 + 2*as(i)*(as(i)*gder*der - 1) / ...
            (as(i)*gval - rlam)^3 * awhts(i);
    end

        end
%
%
%
%
%
        function [val,der] = mpbdry_fminder(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates the minimum of value of F(lambda,e) (for fixed rlam),
%        and the derivative of the minimum as a function of lambda
%
        nsteps = 1000;
        evars = zeros(1,nsteps);
        fders = zeros(1,nsteps);
        fders2 = zeros(1,nsteps);

%
%        initialize to the left of the minimizer (where fder < 0)
%
        emin = mpbdry_leftend(rlam,as,bs,awhts,bwhts,m,n,gam);
        einit=emin+abs(emin) + 10;

        ifder2=0;
        for i=1:1000
%
        [fder,fder2] = mpbdry_evalfder2(einit,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder2);

        if (fder > 0)
%
        einit = (einit + emin)/2;
    end
        if (fder <= 0)
%
        break;
    end
    end

%
%        find the minimum using Newton
%
        evars(1) = einit;
        ifder2=1;
        kstop=0;
        deps = mpbdry_machzero();
        tol=sqrt(deps)/100;
%
        for i=1:nsteps
%
        [fders(i),fders2(i)] = mpbdry_evalfder2(evars(i),rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder2);
        evars(i+1) = evars(i) - fders(i) / fders2(i);
%
        if (abs(fders(i)) < tol)
%
        kstop = kstop+1;
    end
        if (kstop == 2)
%
        nsteps = i;
        break;
    end
    end

        fders=fders(1:nsteps);
        fders2=fders2(1:nsteps);
        evars=evars(1:nsteps);

        val = evars(nsteps);
%
%        evaluate the derivative of the minimum, wrt lambda
%
        [fder2,fdmix] = mpbdry_evalfmix(val,rlam,as,bs,awhts,bwhts,m,n,gam);
        der = -fdmix / fder2;


        end
%
%
%
%
%
        function [fder2,fdmix] = mpbdry_evalfmix(evar,rlam,as,bs,...
            awhts,bwhts,m,n,gam)
%
%        Evaluates the second derivative of F(lambda,e) wrt e, and its 
%        mixed partial derivative
%
%        . . . G and its first derivative (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);
%
%        second derivative of F (wrt e)
%
        term1 = sum(bs.^3 .* bwhts ./ (1+gam*evar*bs).^3);
        term2 = sum(awhts .* as.^2 ./ (as*gval - rlam).^2);
        term3 = sum(awhts .* as.^3 ./ (as*gval - rlam).^3);

        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;
%
%        mixed partial of F
%
        fdmix = 2*gder*sum(awhts .* as.^2 ./ (as*gval - rlam).^3);


        end
%
%
%
%
%
        function [fval,fder,fder2,fder3,gval,gder,gder2,gder3,fdlam,fdmix] = ...
            mpbdry_evalfg3(evar,rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates F(lambda,e), its first three partial derivatives with
%        respect to e; its derivative wrt lambda; its mixed partial; and 
%        G(e) and its first three derivatives
%
%
%        . . . G and its derivatives (wrt e)
%
        gval = sum(bs .* bwhts./ (1 + gam*evar*bs));
        gder = -gam*sum(bs.^2 .* bwhts./ (1 + gam*evar*bs).^2);
        gder2 = 2*gam^2*sum(bs.^3 .* bwhts./ (1 + gam*evar*bs).^3);
        gder3 = -6*gam^3*sum(bs.^4 .* bwhts./ (1 + gam*evar*bs).^4);

%
%        derivative of F wrt lambda
%
        fdlam = -sum(as .* awhts./ (as.*gval - rlam).^2);

%
%        F and its first derivative (wrt e)
%
        fval = evar - sum(as.*awhts ./ (as*gval - rlam));
        fder = 1 + gder*sum(awhts.*as.^2 ./ (as*gval - rlam).^2);

%
%        second and third derivatives of F (wrt e)
%
        term1 = sum(bs.^3 .* bwhts ./ (1+gam*evar*bs).^3);
        dterm1 = -3*sum(bs.^4 .* bwhts ./ (1+gam*evar*bs).^4);

        term2 = sum(awhts .* as.^2 ./ (as*gval - rlam).^2);
        dterm2 = -2*gder*sum(awhts .* as.^3 ./ (as*gval - rlam).^3);


        term3 = sum(awhts .* as.^3 ./ (as*gval - rlam).^3);
        dterm3 = -3*gder*sum(awhts .* as.^4 ./ (as*gval - rlam).^4);


        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;
        fder3 = 2*gam^2*(term1*dterm2 + dterm1*term2) ...
            - 2*(2*gder*gder2*term3 + gder^2*dterm3);

%
%        mixed partial of F
%
        fdmix = 2*gder*sum(awhts .* as.^2 ./ (as*gval - rlam).^3);


        end
%
%
%
%
%
        function [fval,fder,fder2,fder3,gval,gder,gder2,gder3,fdlam,fdmix] = ...
            mpbdry_evalfg3_slow(evar,rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates F(lambda,e), its first three partial derivatives with
%        respect to e; its derivative wrt lambda; its mixed partial; and 
%        G(e) and its first three derivatives
%
%
%        . . . G and its derivatives (wrt e)
%
        gval = 0;
        gder=0;
        gder2 =0;
        gder3=0;
        for j=1:n
%
        gval = gval + bs(j) / (1 + gam*bs(j)*evar) * bwhts(j);
        gder = gder - gam*bs(j)^2 / (1 + gam*bs(j)*evar)^2  * bwhts(j);
        gder2 = gder2 + 2*gam^2*bs(j)^3 / (1 + gam*bs(j)*evar)^3  * bwhts(j);
        gder3 = gder3 - 6*gam^3*bs(j)^4 / (1 + gam*evar*bs(j))^4 * bwhts(j);
    end

%
%        derivative of F wrt lambda
%
        fdlam = 0;
        for j=1:m
%
        fdlam = fdlam - as(j) / (as(j)*gval - rlam)^2 * awhts(j);
    end

%
%        F and its first derivative (wrt e)
%
        fval = 0;
        fder = 0;
        for i=1:m
%
        fval = fval + as(i) / (as(i)*gval - rlam) * awhts(i);
        fder = fder - gder*as(i)^2 / (as(i)*gval - rlam)^2  * awhts(i);
    end
        fval = evar - fval;
        fder = 1-fder;

%
%        second and third derivatives of F (wrt e)
%
        term1 = 0;
        term2 = 0;
        term3 = 0;

        dterm1=0;
        dterm2=0;
        dterm3=0;

        for i=1:n
%
        term1 = term1+bs(i)^3/(1 + gam*bs(i)*evar)^3 * bwhts(i);
        dterm1 = dterm1-3*gam*bs(i)^4/(1 + gam*bs(i)*evar)^4 * bwhts(i);
    end
%
        for i=1:m
%
        term2 = term2 + as(i)^2 / (as(i) * gval - rlam)^2 * awhts(i);
        dterm2 = dterm2 -2*gder* as(i)^3 / (as(i) * gval - rlam)^3 * awhts(i);

        term3 = term3 + as(i)^3 / (as(i)*gval - rlam)^3 * awhts(i);
        dterm3 = dterm3 - 3*gder*as(i)^4 / (as(i)*gval - rlam)^4 * awhts(i);
    end
%
        fder2 = 2*gam^2*term1*term2 - 2*gder^2*term3;
        fder3 = 2*gam^2*(term1*dterm2 + dterm1*term2) ...
            - 2*(2*gder*gder2*term3 + gder^2*dterm3);

%
%        mixed partial of F
%
        fdmix=0;
        for i=1:m
%
        fdmix = fdmix + 2*gder*as(i)^2 / (as(i)*gval - rlam)^3 * awhts(i);
    end


        end
%
