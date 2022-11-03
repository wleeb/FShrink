        function iout = prini(num,init)
        persistent iprin
%
        if (num<0)
%
        iout = iprin;
        return;
    end

        if (init == 1)
%
        prinset(num);
    end

        iprin=num;
        iout=iprin;

        end
%
