        function prod = whtd_scapro(x,y,n)
%
        prod=sum(x(1:n).*y(1:n));

        return;

        x=x(:);
        y=y(:);

        prod = 0;
        for i=1:n
%
        prod = prod + x(i)*y(i);
    end
        end
%
