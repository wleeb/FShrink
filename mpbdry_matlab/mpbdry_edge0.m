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
