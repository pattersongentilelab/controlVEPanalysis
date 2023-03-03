% control subjects VEPs
function [T]=calcVEPpeak(xdata,ydata)


    %% Calculate N75, P100, and N135 peak latency and N75 to P100 peak-to-peak amplitudes and plot results (based on ISCEV 2016 standards)

    
    T=table('Size',[size(ydata,1) 8],'VariableTypes',{'double',...
        'double','double','double','double','double','double','double'},...
        'VariableNames',{'N1_latency','P1_latency','N2_latency',...
        'N1_amplitude','P1_amplitude','N2_amplitude','N1P1_amplitude','P1N2_amplitude'});
    
    for x=1:size(ydata,1)
        temp_ydata=ydata(x,:);
        
        N1=find(temp_ydata==min(temp_ydata(:,62:93)));
        N1=N1(1);

        P1=find(temp_ydata==max(temp_ydata(:,93:133)));
        P1=P1(1);
        
        N2=find(temp_ydata==min(temp_ydata(:,123:179)));
        N2=N2(1);
        
        T(x,:).N1_latency=xdata(:,N1);
        T(x,:).P1_latency=xdata(:,P1);
        T(x,:).N2_latency=xdata(:,N2);
        T(x,:).N1_amplitude=temp_ydata(:,N1);
        T(x,:).P1_amplitude=temp_ydata(:,P1);
        T(x,:).N2_amplitude=temp_ydata(:,N2);
        T(x,:).N1P1_amplitude=abs(diff([temp_ydata(:,N1) temp_ydata(:,P1)]));
        T(x,:).P1N2_amplitude=abs(diff([temp_ydata(:,N2) temp_ydata(:,P1)]));
        
%         figure(150)
%         hold on
%         plot(xdata,temp_ydata,'-k')
%         plot(xdata(:,N1),temp_ydata(:,N1),'+r')
%         plot(xdata(:,P1),temp_ydata(:,P1),'+r')
%         plot(xdata(:,N2),temp_ydata(:,N2),'+r')
%         plot(0.01,abs(diff([temp_ydata(:,N1) temp_ydata(:,P1)])),'+b')
%         plot(0.01,abs(diff([temp_ydata(:,N2) temp_ydata(:,P1)])),'+b')
%         ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0 0.5];ax.YLim=[-1 1];
%         title(num2str(x))
%         pause
%         clf
    end
    
    [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,ydata,...
            'errorbars','Boot','color_mean',[0 0 0],'color_err',[0.8 0.8 0.8],'fig_num',50,...
        'sub_plot',true,'sub_plot_num',[1 2 1]);
    
    hold on
    plot(T.N1_latency,-0.11,'ok')
    plot(T.P1_latency,0.04,'ok')
    plot(T.N2_latency,-0.05,'ok')
    ax=gca;ax.YLim=[-0.15 0.15];
    
    subplot(1,2,2)
    hold on
    x1=mean([mean(T.N1_latency) mean(T.P1_latency)]);
    x2=mean([mean(T.N2_latency) mean(T.P1_latency)]);
    plot(x1,T.N1P1_amplitude,'ok')
    plot(x2,T.P1N2_amplitude,'ok')
    ax=gca;ax.Box='off';ax.TickDir='out';ax.XTick=[x1 x2];...
    ax.XTickLabel={'N1P1','P1N2'};ax.XLim=[0.05 0.15];ax.YLim=[0 0.3];
    
end