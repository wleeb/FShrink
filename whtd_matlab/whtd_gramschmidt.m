        function u = whtd_gramschmidt(u,m,k)
%
        u(:,1) = u(:,1) / norm(u(:,1));

        for i=2:k

        for ijk=1:2
%
        for j=1:i-1
%
        pij = sum(u(:,i) .* u(:,j));
        u(:,i) = u(:,i) - pij*u(:,j);
    end
        u(:,i) = u(:,i) / norm(u(:,i));
    end
    end

        end
%
