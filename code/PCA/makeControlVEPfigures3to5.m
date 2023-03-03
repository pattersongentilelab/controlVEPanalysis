% Fit control subject data to PCA model

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/cleaned_VEP.mat'])

xdata=cell2mat(cleaned_vep(1,3));
PC_no=7;
edges=-2.5:0.25:2.5;
center=edges(1:end-1)+diff(edges);
X=-2.5:0.01:2.5;

%% Subject selection

% select control subjects
control_subject=find(cleaned_vep_files.subjecttype=='Control');
cleaned_vep=cleaned_vep(control_subject,:);
control_vep_subject_sessions=cleaned_vep_files(control_subject,:);

unique_ID=unique(control_vep_subject_sessions.uniqueID);

%% Load PCA model
data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/PCAmodel.mat'])

%% Fit PCA model to the whole dataset

for x=1:length(unique_ID)
    temp_loc=find(cell2mat(cleaned_vep(:,1))==unique_ID(x,:));
    if length(temp_loc)>1
        temp_loc=temp_loc(1:2);
        mAge(x,:)=mean(control_vep_subject_sessions.age_vep(temp_loc,:));
    else
        mAge(x,:)=control_vep_subject_sessions.age_vep(temp_loc,:);
        
    end
    temp_ydata=cleaned_vep(temp_loc,4);
    control_vep_subjects(x,:)=control_vep_subject_sessions(temp_loc(1),:);
    for y=1:size(temp_ydata,1)
        temp_ydata2(y,:)=nanmean(cell2mat(temp_ydata(y,:)));
    end
    switch size(temp_ydata2,1)
        case 1
            controlM_vep(x,:)=temp_ydata2;
        otherwise
            controlM_vep(x,:)=nanmean(temp_ydata2,1);
    end
    clear temp_ydata*

end

control_vep_subjects.age_vep=mAge;

for x=1:size(controlM_vep,1)
    for y=1:size(train_coeff,2)
        temp1=controlM_vep(x,:);
        temp2=train_coeff(:,y);
        r=temp1*temp2;
        controlM_score(x,y)=r;
    end
end


  %% ANOVA calculation on subgroups 
  % subgroups med Hx:concussion history, family history of migraine, history of migraine
  % male vs. female
  migraine_Hx=control_vep_subjects.med_hx___1;
  migraine_famHx=control_vep_subjects.family_hx___1;
  concussion_Hx=control_vep_subjects.conc_hx_yn;
  sex=control_vep_subjects.sex_master;
  age=control_vep_subjects.age_vep;

  %% Adjust PC scores so all subjects are sex nulled and nulled to age 15 years
[ADJcontrolM_score]=AgeSexNull(7,controlM_score,sex,age);

  %% Adjust PC scores so all subjects are nulled to age 15 years
[ADJAcontrolM_score]=AgeSexNull(7,controlM_score,[],age);

  %% Adjust PC scores so all subjects are sex nulled
[ADJScontrolM_score]=AgeSexNull(7,controlM_score,sex,[]);

%% age plot
figure(10)
% plot PC2
subplot(1,2,1)
hold on
plot(control_vep_subjects.age_vep,ADJScontrolM_score(:,2),'.','Color',[0.5 0.5 0.5],'MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[-2.2 2.2];
axis('square')
[corr p_val]=corrcoef(control_vep_subjects.age_vep,ADJScontrolM_score(:,2));
title(['PC2, corrcoef=' num2str(corr(1,2)) ', p-val=' num2str(p_val(1,2))])
xlabel('Age at time of first VEP')
ylabel('Component score')


% plot PC4
subplot(1,2,2)
hold on
plot(control_vep_subjects.age_vep,ADJScontrolM_score(:,4),'.','Color',[0.5 0.5 0.5],'MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[-2.2 2.2];
axis('square')
[corr p_val]=corrcoef(control_vep_subjects.age_vep,ADJScontrolM_score(:,4));
title(['PC4, corrcoef=' num2str(corr(1,2)) ', p-val=' num2str(p_val(1,2))])
xlabel('mean age across VEP sessions')
ylabel('Component score')


%% compare sexes
% female_subjects=find(control_vep_subjects.sex_master==1 & control_vep_subjects.med_hx___1~=1);
% male_subjects=find(control_vep_subjects.sex_master==2 & control_vep_subjects.med_hx___1~=1);

female_subjects=find(control_vep_subjects.sex_master==1);
male_subjects=find(control_vep_subjects.sex_master==2);


 [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,controlM_vep(female_subjects,:),...
                'errorbars','Boot','color_mean',[0 0 0],'color_err',[0.5 0.5 0.5],'fig_num',11,...
                'sub_plot',true,'sub_plot_num',[1 3 1]);
        
  [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,controlM_vep(male_subjects,:),...
                'errorbars','Boot','color_mean',[0 1 1],'color_err',[0.5 1 1],'fig_num',11,...
                'sub_plot',true,'sub_plot_num',[1 3 1]);
            
  legend('show')
  legend('females','','males','')

% PC1
subplot(1,3,2)
hold on
hist1=histcounts(ADJAcontrolM_score(female_subjects,1),edges,'Normalization','probability');
Y1=spline(center,hist1,X);
hist2=histcounts(ADJAcontrolM_score(male_subjects,1),edges,'Normalization','probability');
Y2=spline(center,hist2,X);
plot(center,hist1,'.k',X,Y1,'-k')
plot(center,hist2,'.c',X,Y2,'-c')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-2.5 2.5];ax.YLim=[0 0.5];
plot(nanmedian(ADJAcontrolM_score(female_subjects,1),1),0.45,'vk')
plot(nanmedian(ADJAcontrolM_score(male_subjects,1),1),0.45,'vc')
[h p]=ttest2(ADJAcontrolM_score(female_subjects,1),ADJAcontrolM_score(male_subjects,1));
title(['PC1, p-value=' num2str(p)])
xlabel('Component score')
ylabel('Number of subjects')
legend('show')
legend('Female','Male')

% PC4
subplot(1,3,3)
hold on
hist1=histcounts(ADJAcontrolM_score(female_subjects,4),edges,'Normalization','probability');
Y1=spline(center,hist1,X);
hist2=histcounts(ADJAcontrolM_score(male_subjects,4),edges,'Normalization','probability');
Y2=spline(center,hist2,X);
plot(center,hist1,'.k',X,Y1,'-k')
plot(center,hist2,'.c',X,Y2,'-c')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-2.5 2.5];ax.YLim=[0 0.5];
plot(nanmedian(ADJAcontrolM_score(female_subjects,4),1),0.45,'vk')
plot(nanmedian(ADJAcontrolM_score(male_subjects,4),1),0.45,'vc')
[h p]=ttest2(ADJAcontrolM_score(female_subjects,4),ADJAcontrolM_score(male_subjects,4));
title(['PC4, p-value=' num2str(p)])
xlabel('Component score')
ylabel('Number of subjects')
legend('show')
legend('Female','Male')

figure(20)
subplot(1,2,1)
hold on
for x=1:PC_no
  [dprime_sex(x,:) pval_sex(x,:)]=calcDprime(ADJAcontrolM_score(female_subjects,x),ADJAcontrolM_score(male_subjects,x));
end
errorbar(dprime_sex(:,2),1:PC_no,[],[],abs(diff(dprime_sex(:,1:2),[],2)),abs(diff(dprime_sex(:,2:3),[],2)),'ok')
ax=gca;ax.Box='off';ax.TickDir='out';ax.YLim=[0 PC_no+1];ax.XLim=[-1 1.1];ax.YDir='reverse';
plot([0 0],[0 PC_no+1],'--k')
xlabel('Effect size')
ylabel('Component')
title('sex')

%% Simulated VEP for a 10 year old, 15 year old, and 20 year old sex-nulled subject

[simVEP]=synthVEP(7,ADJScontrolM_score,train_coeff,explained','age',age);
