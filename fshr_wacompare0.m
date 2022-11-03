        function fshr_wacompare0(evals,gam,gam0,ifig)
%
        [x0,x1] = fshr_wachter_lims(gam,gam0);

        nvals=500;
        vals = zeros(1,nvals);
        ts=linspace(x0,x1,nvals);

        for i=1:nvals
%
        vals(i) = fshr_evalwacht(ts(i),gam,gam0);
    end

        figure(ifig);
        subplot(1,2,1);
        hh = histogram(evals,100,'Normalization','pdf');
        hold on; plot(ts,vals,'linewidth',3)

        m=length(evals);

        subplot(1,2,2);
        plot(evals,'*','Color','b')
        hold on;
        plot(0:m+1,x1*ones(m+2,1),'LineWidth',2,'Color','r')
        plot(0:m+1,x0*ones(m+2,1),'LineWidth',2,'Color','r')

        xlim([0,m+1])

        set(figure(ifig),'Position',[500,500,1300,500])

        end
