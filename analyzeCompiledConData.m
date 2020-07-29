% analyze the compiled data

filenameMAT=fullfile('/Users','carlynpattersongentile','Documents','MATLAB','concussion VEP','compiled_ConVEP_data.mat');

load(filenameMAT);

clear ans filenameMAT

% plot mean VEP data
for XX=1:size(vep_by_subjectHC,2)
    y_data_subjectsHCm(XX,:)=mean(cell2mat(vep_by_subjectHC(:,XX)),1);
    figure
    
%     figure(10)
%     Y_data=cell2mat(vep_by_subjectHC(:,XX));
%     for z=1:size(Y_data,1)
%         plot(x_data,Y_data(z,:),'b')
%         ax=gca; ax.TickDir='out'; ax.Box='off';
%         title(num2str(XX))
%         ylabel('VEP signal (mV)')
%         xlabel('Time (ms)')
%         pause
%         hold off
%     end
end

plot_meanVEP(x_data,y_data_subjectsHCm,'plot_all',true,'normalize',true,'errorbars','Boot','fig_num',1);



% Look at habituation
for XX=1:size(vep_by_subjectHC,2)
    temp=cell2mat(vep_by_subjectHC(:,XX));
    y_data_subjectsHCm_early(XX,:)=mean(temp(1:8,:),1);
end

plot_meanVEP(x_data,y_data_subjectsHCm_early,'plot_mean',true,'normalize',true,'errorbars','Boot','fig_num',2);

for XX=1:size(vep_by_subjectHC,2)
    temp=cell2mat(vep_by_subjectHC(:,XX));
    y_data_subjectsHCm_late(XX,:)=mean(temp(end-8:end,:),1);
end

plot_meanVEP(x_data,y_data_subjectsHCm_late,'plot_mean',true,'normalize',true,'errorbars','Boot','fig_num',2,'color_mean',[1 0 0],'color_err',[1 0.8 0.8]);



% extract gamma oscillations
for XX=1:size(vep_by_subjectHC,2)
    disp(num2str(XX))
    y_data=cell2mat(vep_by_subjectHC(:,XX));
    [cfs,f,psd,freqs]=calc_plot_TimeFreqDomain(max(x_data)/2*1000,Fs,x_data,y_data,'normalize',true,'induced',true,'plot_all',true);
    CFS(XX)={cfs};
    CFSm(XX,:,:)=squeeze(mean(cfs,1));
    clear y_data
end

% gamma_temp=[2:6 12 15 24 41 45 47 68 70 71 87 96 120 131];
% gamma_temp=[5 9 13 14 17 20 21 25 27 29];

figure
cfsM=squeeze(mean(CFSm,1));
imagesc(x_data,log2(f),cfsM);
Yticks=2.^(round(log2(min(f))):round(log2(max(f))));
ax=gca;
ax.YLim=[min(f) max(f)];
ax.YTick=log2(Yticks);
ax.YDir='normal';
set(ax,'YLim',log2([min(f),max(f)]), ...
'layer','top', ...
'YTick',log2(Yticks(:)), ...
'YTickLabel',num2str(sprintf('%g\n',Yticks)), ...
'layer','top')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
title('CWT of prVEP signal');
hcol=colorbar;
cmap=colormap('jet');
cmap(29:36,:)=0.5;
colormap(cmap)
caxis([-2 2])
hcol.Label.String='Normalized Magnitude (z-score)';
