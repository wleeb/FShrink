        function main
        prini(13,1)


        n=16
        x = randn(1,n);
        y = randn(1,n);
        z = randn(1,n);


        prin2('x=',x,n);


        prini(0,0);
        prin2('y=',y,n)

        aaaa=prini(13,0);
        prin2('z=',z,n)




        mms=randi(200,1,20);

        mms = mms.*sign(randn(1,20));


        prinf('mms=',mms,20);

        a = randn(7,12)


        prinar2('a=',a,4,2)

        prinar2('a=',a,5,6)




        mar = randi(200000,12,15).*sign(randn(12,15));

        prinarf('mar=',mar,6,10)

        prinstop
        end
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This code contains several user-callable functions for printing.
%
%    prini - initializes the print routine and turns printing on/off
%    prin2 - prints 1-dimensional arrays of doubles
%    prinar2 - prints 2-dimensional arrays of doubles
%    prinf - prints 1-dimensional integer arrays
%    prinarf - prints 2-dimensional integer arrays
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
        function prinarf(msg,x,m,n)
%
        iout = prini(-1,0);
        if (iout <= 0)
%
        return
    end
%
        x=x(1:m,1:n)';

        fmt = repmat('  %7d',1,n);
        fmt=strcat(['  \n',fmt]);

        fprintf(strcat(['  ',msg]));
        fprintf(fmt,x);
        fprintf('\n\n');


        end
%
%
%
%
%
        function prinar2(msg,x,m,n)
%
        iout = prini(-1,0);
        if (iout <= 0)
%
        return
    end
%
        x=x(1:m,1:n)';

        fmt = repmat('    %+3.5e',1,n);
        fmt=strcat(['  \n',fmt]);

        fprintf(strcat(['  ',msg]));
        fprintf(fmt,x);
        fprintf('\n\n');


        end
%
%
%
%
%
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
%
%
%
%
        function prin2(msg,x,n)
%
        iout = prini(-1,0);
        if (iout <= 0)
%
        return
    end
%
        len = 6;
        mmm = mod(n,len);
        m = (n - mmm) / len;
%
        fmt0 = repmat('    %+3.5e',1,len);
        fmt0 = repmat([fmt0,'\n'],1,m);
        fmt1 = repmat('    %+3.5e',1,mmm);
        fmt = strcat(['  ',msg,'\n',fmt0,fmt1,'\n\n']);

        fprintf(fmt,x(1:n));

        end
%
%
%
%
%
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
%
%
%
%
        function prinset(num)
%
        fname = strcat('out',num2str(num));
        delete(fname)
        diary(fname)
        diary on
%
        format short E
        rng('default');

        end
%
%
%
%
%
        function prinstop()
%
        diary off
        error('stop')

        end
%
%
