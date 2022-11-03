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
