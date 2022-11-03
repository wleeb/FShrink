        function a2 = whtd_dmult(a,dl,dr,m,n)
%
%        forms the product
%
%             A2 = Dleft * A * Dright,                        (1)
%
%        where Dleft and Dright are specified diagonal matrices
%
%
%        . . . ensure dl and dr are column vectors
%
        [n1,n2] = size(dl);
        if (n1 < n2)
%
        dl=dl';
    end
        [n1,n2] = size(dr);
        if (n1 < n2)
%
        dr=dr';
    end
%
%        form the product
%
        a2 = bsxfun(@times,a,dl);
        a2 = bsxfun(@times,a2',dr)';

        end
%
