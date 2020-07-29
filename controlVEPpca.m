% control subjects VEPs

load cleaned_VEP control_vep

% organize control VEP ydata

for x=1:size(control_vep,1)
    temp_ydata=cell2mat(control_vep(x,2));
    control_ydata(x,:)=nanmedian(temp_ydata,1);
end

xdata=cell2mat(control_vep(1,3));

[coeff,score,latent,tsquared,explained,mu]=pca(control_ydata);


[y_dataM,y_dataERR1,y_dataERR2]=plot_medianVEP(xdata,control_ydata,...
            'errorbars','Boot','color_mean',[0 0 0],'color_err',[0.8 0.8 0.8],'fig_num',1,...
        'sub_plot',true,'sub_plot_num',[1 2 1]);

subplot(1,2,2)
hold on
plot(xdata,coeff(:,1:3)*-1)
legend('show')
legend('PC1','PC2','PC3','PC4','PC5')
title([ 'Variance explained (%) = ' num2str(sum(explained(1:3)))])
ax=gca;ax.Box='off';ax.TickDir='out';