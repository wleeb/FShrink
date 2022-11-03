        function a2 = fshr_dleft(d,a,m,n)
%
%        forms the product
%
%             A2 = D * A,                               (1)
%
%        where D is a specified diagonal matrix
%

%
%        ensure diag is a column vector
%
        [n1,n2] = size(d);
        if (n1 < n2)
%
        d=d';
    end
%
%        form the product
%
        a2 = bsxfun(@times,a,d);

        end
