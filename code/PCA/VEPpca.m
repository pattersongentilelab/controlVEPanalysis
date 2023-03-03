% control subjects VEPs
function [coeff,score,explained,latent]=VEPpca(xdata,ydata,PC_no)


    %% Calculate PCA and plot results

    % get one averaged waveform per subject across all sessions, case and
    % control subjects

    [coeff,score,latent,tsquared,explained,mu]=pca(ydata,'Centered',false);

    [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,ydata,...
                'errorbars','Boot','color_mean',[0 0 0],'color_err',[0.8 0.8 0.8],'fig_num',1,...
            'sub_plot',true,'sub_plot_num',[3 1 1]);
        
   % Determine PC sign
   leg_label={};
   for x=1:size(coeff,2)
       temp=nanmedian(score(:,x),1);
       leg_label=cat(2,leg_label,{['PC' num2str(x)]});
       if temp >=0
           sign_pc(:,x)=1;
       else
           sign_pc(:,x)=-1;
       end
   end
   sign_pc1=repmat(sign_pc,length(xdata),1);
   coeff=coeff.*sign_pc1;
   sign_pc2=repmat(sign_pc,size(score,1),1);
   score=score.*sign_pc2;

    subplot(3,1,2)
    hold on
    for i=1:PC_no
        plot(xdata,coeff(:,i).*latent(i,:))
    end
    legend('show')
    legend(leg_label)
    title([ 'Variance explained (%) = ' num2str(sum(explained(1:PC_no)))])
    ax=gca;ax.Box='off';ax.TickDir='out';

    subplot(3,1,3)
    pareto(explained)
    ax=gca;ax.Box='off';ax.TickDir='out';

    % Plot scores of PCs as histograms
    figure(2)

    for x=1:PC_no
        subplot(ceil(PC_no/3),3,x)
        histogram(score(:,x),-1.6:0.2:1.6)
        ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-1.6 1.6];
        title(['Principal component ' num2str(x)])
        if x==4
            xlabel('Component score')
        end
        if x==PC_no
            ylabel('Number of subjects')
        end
    end


end