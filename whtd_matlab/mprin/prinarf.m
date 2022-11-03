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
