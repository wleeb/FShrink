        function prinf(msg,x,n)
%
        iout = prini(-1,0);
        if (iout <= 0)
%
        return
    end
%
        len = 10;
        mmm = mod(n,len);
        m = (n - mmm) / len;
%
        fmt0 = repmat('    %5d',1,len);
        fmt0 = repmat([fmt0,'\n'],1,m);
        fmt1 = repmat('    %5d',1,mmm);
        fmt = strcat(['  ',msg,'\n',fmt0,fmt1,'\n\n']);

        fprintf(fmt,x(1:n));

        end
%
