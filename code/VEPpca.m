% control subjects VEPs

load cleaned_VEP

%% Subject selection

% select only control subjects
control_subject=find(cleaned_vep_files.subjecttype=='Control');
cleaned_vep=cleaned_vep(control_subject,:);
cleaned_vep_files=cleaned_vep_files(control_subject,:);

unique_ID=unique(cleaned_vep_files.uniqueID);


counter=1;
for x=1:length(unique_ID)
    temp_loc=find(cell2mat(cleaned_vep(:,1))==unique_ID(x,:));
    if length(temp_loc)>1 % select only subjects with multiple sessions
        temp_loc=temp_loc(1:2); % select only the first two sessions
        temp_ydata=cleaned_vep(temp_loc,4);
        for y=1:size(temp_ydata,1)
            cleaned_vepM(counter,:)=nanmean(cell2mat(temp_ydata(y,:)));
            subject_loc(counter,:)=temp_loc(1);
            subject_session_loc(counter,:)=temp_loc(y);
            subject_session_no(counter,:)=length(temp_loc);
            counter=counter+1;
        end
    end
end

subject_loc=unique(subject_loc);

% select only the first 50 subjects to save the rest for validations
cleaned_vepM=cleaned_vepM(1:100,:);
subject_loc=subject_loc(1:50,:);
subject_session_loc=subject_session_loc(1:100,:);
cleaned_vep_subject_files_forM=cleaned_vep_files(subject_loc,:);
cleaned_vep_session_files_forM=cleaned_vep_files(subject_session_loc,:);

%% Calculate PCA and plot results

% get one averaged waveform per subject across all sessions, case and
% control subjects

xdata=cell2mat(cleaned_vep(1,3));

[coeff,score,latent,tsquared,explained,mu]=pca(cleaned_vepM,'Centered',false);


[y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,cleaned_vepM,...
            'errorbars','Boot','color_mean',[0 0 0],'color_err',[0.8 0.8 0.8],'fig_num',1,...
        'sub_plot',true,'sub_plot_num',[3 1 1]);


sign_pc=repmat([-1 1 1 -1 1 1 -1 -1],512,1);    
subplot(3,1,2)
hold on
plot(xdata,coeff(:,1:8).*sign_pc)
legend('show')
legend('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8')
title([ 'Variance explained (%) = ' num2str(sum(explained(1:8)))])
ax=gca;ax.Box='off';ax.TickDir='out';

subplot(3,1,3)
pareto(explained)
ax=gca;ax.Box='off';ax.TickDir='out';

% Plot scores of PCs as histograms
figure(2)

for x=1:8
    subplot(4,2,x)
    histogram(score(:,x),-1.6:0.2:1.6)
    ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-1.6 1.6];
    title(['Principal component ' num2str(x)])
    if x==5 || x==6
        xlabel('Component score')
    end
    if x==3
        ylabel('Number of subjects')
    end
end

%% Determine intra-subject correlation coefficients for each PC
score_session1=score(1:2:end,:);
score_session2=score(2:2:end,:);

figure(3)
for x=1:size(score,2)
    temp_corr=corrcoef(score_session1(:,x),score_session2(:,x));
    corr(x,:)=temp_corr(1,2);
end
plot(corr,'xb')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[1 100];ax.YLim=[-1 1];
axis('square')
xlabel('Principal component')
ylabel('Pearson correlation coefficient (session 1 vs. session 2)')

figure(4)
for x=1:8
    subplot(4,2,x)
    plot(score_session1(:,x),score_session2(:,x),'xb')
    hold on
    plot([-1.6 1.6],[-1.6 1.6],'--k')
    ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-1.6 1.6];ax.YLim=[-1.6 1.6];
    axis('square')
    title(['PC' num2str(x) ' corrcoef=' num2str(corr(x,:))])
    if x==7 
        xlabel('Component score session 1')
    end
    if x==3
        ylabel('Component score session 2')
    end
end


