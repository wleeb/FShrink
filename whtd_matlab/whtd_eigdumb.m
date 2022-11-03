        function [ua,sa] = whtd_eigdumb(a,m)
%
        [ua,sa]=eig(a);
        sa=diag(sa);
        [sa,ii] = sort(sa,'descend');
        ua=ua(:,ii);

        end
%
