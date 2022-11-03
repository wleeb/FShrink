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
