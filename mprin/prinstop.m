        function prinstop()
%
        diary off
%%%        error('stop')

%%%        return;
%
%        safely stop without error message
%
        ms.message='';
        ms.stack = dbstack('-completenames');
        ms.stack(1:end) = [];
        ds = dbstatus();
        stoponerror = any(strcmp('error', {ds.cond}));
        setappdata(0, 'dberrorkeep', stoponerror);
        dbclear error

        error(ms);

        end
%
