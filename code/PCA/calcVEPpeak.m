% control subjects VEPs
function [T]=calcVEPpeak(xdata,ydata)


    %% Calculate N75, P100, and N135 peak latency and N75 to P100 peak-to-peak amplitudes and plot results (based on ISCEV 2016 standards)

    
    T=table('Size',[size(ydata,1) 5],'VariableTypes',{'double',...
        'double','double','double','double'},'VariableNames',{'N1_latency',...
        'P1_latency','N2_latency','N1P1_amplitude','P1N2_amplitude'});
    
    for x=1:size(ydata,1)
        temp_ydata=ydata(x,:);
        N1=find(temp_ydata==min(temp_ydata(:,51:92)));
        if length(N1)>1
            temp=find(N1>50 & N1<=93);
            N1=N1(temp(1));
        end
        P1=find(temp_ydata==max(temp_ydata(:,93:123)));
        if length(P1)>1
            temp=find(P1>92 & P1<=123);
            P1=P1(temp(1));
        end
        N2=find(temp_ydata==min(temp_ydata(:,124:185)));
        if length(N2)>1
            temp=find(N2>123 & N2<=185);
            N2=N2(temp(1));
        end
        T(x,:).N1_latency=xdata(:,N1);
        T(x,:).P1_latency=xdata(:,P1);
        T(x,:).N2_latency=xdata(:,N2);
        T(x,:).N1P1_amplitude=abs(diff([temp_ydata(:,N1) temp_ydata(:,P1)]));
        T(x,:).P1N2_amplitude=abs(diff([temp_ydata(:,N2) temp_ydata(:,P1)]));
    end
    
    [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,ydata,...
            'errorbars','Boot','color_mean',[0 0 0],'color_err',[0.8 0.8 0.8],'fig_num',5,...
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
    
    figure(6)
    subplot(3,2,1)
    temp=corrcoef(T(1:2:end,:).N1_latency,T(2:2:end,:).N1_latency);
    plot(T(1:2:end,:).N1_latency,T(2:2:end,:).N1_latency,'ok')
    hold on
    plot([0.045 0.1],[0.045 0.1],'--k')
    axis('square')
    title(['N1 latency, corr=' num2str(temp(1,2))])
    ax=gca;ax.Box='off';ax.TickDir='out';
    
    subplot(3,2,3)
    temp=corrcoef(T(1:2:end,:).P1_latency,T(2:2:end,:).P1_latency);
    plot(T(1:2:end,:).P1_latency,T(2:2:end,:).P1_latency,'ok')
    hold on
    plot([0.08 0.12],[0.08 0.12],'--k')
    axis('square')
    title(['P1 latency, corr=' num2str(temp(1,2))])
    ax=gca;ax.Box='off';ax.TickDir='out';
    
    subplot(3,2,5)
    temp=corrcoef(T(1:2:end,:).N2_latency,T(2:2:end,:).N2_latency);
    plot(T(1:2:end,:).N2_latency,T(2:2:end,:).N2_latency,'ok')
    hold on
    plot([0.1 0.2],[0.1 0.2],'--k')
    axis('square')
    title(['N2 latency, corr=' num2str(temp(1,2))])
    ax=gca;ax.Box='off';ax.TickDir='out';
    xlabel('session 1')
    
    subplot(3,2,2)
    temp=corrcoef(T(1:2:end,:).N1P1_amplitude,T(2:2:end,:).N1P1_amplitude);
    plot(T(1:2:end,:).N1P1_amplitude,T(2:2:end,:).N1P1_amplitude,'ok')
    hold on
    plot([0 0.4],[0 0.4],'--k')
    axis('square')
    title(['N1P1 amplitude, corr=' num2str(temp(1,2))])
    ax=gca;ax.Box='off';ax.TickDir='out';
    
    subplot(3,2,4)
    temp=corrcoef(T(1:2:end,:).P1N2_amplitude,T(2:2:end,:).P1N2_amplitude);
    plot(T(1:2:end,:).P1N2_amplitude,T(2:2:end,:).P1N2_amplitude,'ok')
    hold on
    plot([0 0.4],[0 0.4],'--k')
    axis('square')
    title(['P1N2 amplitude, corr=' num2str(temp(1,2))])
    ax=gca;ax.Box='off';ax.TickDir='out';
    ylabel('session 2')
    
end