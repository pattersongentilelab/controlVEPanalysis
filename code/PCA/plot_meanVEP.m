function [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(x_data,y_data,varargin)

q = inputParser;
q.addParameter('normalize',false,@islogical);
q.addParameter('fig_num',1,@isnumeric);
q.addParameter('sub_plot',false,@islogical);
q.addParameter('sub_plot_num',[1 1 1],@isnumeric);
q.addParameter('color_mean',[0 0 0],@isnumeric);
q.addParameter('color_err',[0.8 0.8 0.8],@isnumeric);
q.addParameter('errorbars','STEr',@ischar); % errorbar method
% 'STE' is for standard error of the mean errorbars, 'STD' is standard
% deviation, 'Boot' is bootstrapped at 95% confidence intervals

q.parse(varargin{:});

if q.Results.normalize==1
    for x=1:size(y_data,1)
        temp=max(abs(y_data(x,:)));
        y_data(x,:)=y_data(x,:)./temp;
    end
end

if q.Results.errorbars=='STEr'
    y_dataM=mean(y_data,1);
    y_dataERR1=y_dataM-(std(y_data,[],1)./sqrt(size(y_data,1)-1));
    y_dataERR2=y_dataM+(std(y_data,[],1)./sqrt(size(y_data,1)-1));
end

if q.Results.errorbars=='STDe'
    y_dataM=mean(y_data,1);
    y_dataERR1=y_dataM-(std(y_data,[],1));
    y_dataERR2=y_dataM+(std(y_data,[],1));
end

if q.Results.errorbars=='Boot'
    bootval=bootstrp(1000,@nanmean,y_data);
    bootval=sort(bootval);
    y_dataM=bootval(500,:);
    y_dataERR1=bootval(25,:);
    y_dataERR2=bootval(975,:);
end

x_ERR=cat(2,x_data,fliplr(x_data));
y_ERR=cat(2,y_dataERR1,fliplr(y_dataERR2));

figure(q.Results.fig_num)
if q.Results.sub_plot==1
    subplot(q.Results.sub_plot_num(1),q.Results.sub_plot_num(2),q.Results.sub_plot_num(3))
end
hold on
TEMP=fill(x_ERR,y_ERR,q.Results.color_err,'EdgeColor','none');
plot(x_data,y_dataM,'-','Color',q.Results.color_mean)
ax=gca; ax.TickDir='out'; ax.Box='off';ax.XLim=[0 0.5];
title('Mean VEP signal')
ylabel('VEP signal (microV)')
xlabel('Time (s)')
    if q.Results.normalize==1
        ylabel('Normalized VEP signal')
        ax.YLim=[-1 1];
    end
end
    