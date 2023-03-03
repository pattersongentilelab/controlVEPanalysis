%% 
% control subjects VEPs

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
analysis_path = getpref('controlVEPanalysis','controlVEP_AnalysisPath');
load([data_path '/cleaned_VEP.mat'])

xdata=cell2mat(cleaned_vep(1,3));
%% Subject selection

% select control subjects
control_subject=find(cleaned_vep_files.subjecttype=='Control');
cleaned_vep=cleaned_vep(control_subject,:);
control_vep_subjects=cleaned_vep_files(control_subject,:);

unique_ID=unique(control_vep_subjects.uniqueID);

counter=1;
counter2=1;
for x=1:length(unique_ID)
    temp_loc=find(cell2mat(cleaned_vep(:,1))==unique_ID(x,:));
    subject_loc(x,:)=temp_loc(1);
    if length(temp_loc)>1 % select only subjects with multiple sessions <6 months apart
        age_range=control_vep_subjects.age_vep(temp_loc(1:2),:);
        if diff(age_range,[],1)<0.5
        temp_loc=temp_loc(1:2); % select only the first two sessions
        temp_ydata=cleaned_vep(temp_loc,4);
            for y=1:size(temp_ydata,1)
                control_vep(counter2,y,:)=nanmean(cell2mat(temp_ydata(y,:)));
                counter=counter+1;
            end
        mAge(counter2,:)=mean(control_vep_subjects.age_vep(temp_loc(1:2),:));
        Age_range(counter2,:)=control_vep_subjects.age_vep(temp_loc(1:2),:);
        control_vep_subjects_sess2(counter2,:)=control_vep_subjects(temp_loc(1),:);
        control_vep_subjects_sess22(counter2,:)=control_vep_subjects(temp_loc(2),:);
        counter2=counter2+1;
        end
    end
    clear temp_loc temp_ydata
end

control_vep_subjects_sess2.age_vep=mAge;
age_diff=diff(Age_range,[],2);

figure(200)
hold on
plot(ones(size(age_diff)),age_diff,'.k')

control_vep_subjects=control_vep_subjects(subject_loc,:);

% load randomly selected 20 male and 20 female subjects for the training PCA
female=find(control_vep_subjects_sess2.sex_master==1);
male=find(control_vep_subjects_sess2.sex_master==2);

load([data_path '/randControlTrainTest.mat'])

tempF1=sort(tempF(1:20));
tempF2=sort(tempF(21:end));
tempM1=sort(tempM(1:20));
tempM2=sort(tempM(21:end));

FM1=cat(1,female(tempF1),male(tempM1));
FM2=cat(1,female(tempF2),male(tempM2));

control_train=control_vep_subjects_sess2(FM1,:);
control_train_vep=control_vep(FM1,:,:);
control_pool=control_vep_subjects_sess2(FM2,:);
control_pool_vep=control_vep(FM2,:,:);

age_diff_train=age_diff(FM1,:);


% find age and sex-matched control subject for every control train subject
counter_control_subjects=1:size(control_pool,1);
for x=1:size(control_train,1)
    temp_age=control_train.age_vep(x,:);
    temp_sex=control_train.sex_master(x,:);
    control_age=abs(control_pool.age_vep-temp_age);
    sex_match=find(control_pool.sex_master==temp_sex);
    age_match=find(control_age(sex_match)==min(control_age(sex_match)));
    sub_match=sex_match(age_match);
    sub_match=sub_match(1);
    control_test(x,:)=control_pool(sub_match,:);
    control_test_vep(x,:,:)=control_pool_vep(sub_match,:,:);
    age_diff_test(x,:)=age_diff(sub_match,:);
    counter_control_subjects=1:size(control_pool,1);
    counter_control_subjects=counter_control_subjects(find(counter_control_subjects~=sub_match));
    control_pool=control_pool(counter_control_subjects,:);
    control_pool_vep=control_pool_vep(counter_control_subjects,:,:);
end

save([analysis_path '/controlTrainTest'],'control_train_vep','control_test_vep','control_train','control_test','xdata');

ydata_train=squeeze(nanmean(control_train_vep,2));
ydata_test=squeeze(nanmean(control_test_vep,2));

%% Calculate PCA and plot results on train group
PC_no=7; % principal component number
[train_coeff,score,explained,latent]=VEPpca(xdata,ydata_train,PC_no);

save([analysis_path '/PCAmodel'],'train_coeff','score','explained');

%% fit train PCA to the remaining test control subjects (sex- and age-matched)

for x=1:size(ydata_train,1)
    for y=1:size(train_coeff,2)
        temp1=ydata_train(x,:);
        temp2=train_coeff(:,y);
        r=temp1*temp2;
        train_score(x,y)=r;
    end
end

for x=1:size(ydata_test,1)
    for y=1:size(train_coeff,2)
        temp1=ydata_test(x,:);
        temp2=train_coeff(:,y);
        r=temp1*temp2;
        test_score(x,y)=r;
    end
end

% plot scores
figure(3)
edges=-2.5:0.25:2.5;
center=edges(1:end-1)+diff(edges);
X=-2.5:0.01:2.5;
for x=1:PC_no
    subplot(ceil(PC_no/3),3,x)
    hold on
    hist1=histcounts(train_score(:,x),edges,'Normalization','probability');
    Y1=spline(center,hist1,X);
    hist2=histcounts(test_score(:,x),edges,'Normalization','probability');
    Y2=spline(center,hist2,X);
    plot(center,hist1,'.k',X,Y1,'-k')
    plot(center,hist2,'.b',X,Y2,'-b')
    [h p]=kstest2(train_score(:,x),test_score(:,x));
    ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-2.5 2.5];ax.YLim=[0 0.5];
    title(['p-value=' num2str(p)]);
end

mdl_test=fitlm(cat(1,train_score(:,PC_no),test_score(:,PC_no)),cat(1,zeros(40,1),ones(40,1)));

anova(mdl_test,'summary')


%% Determine intra-subject correlation coefficients for each PC in train and
% test groups

for x=1:size(control_train_vep,1)
    for y=1:size(control_train_vep,2)
        for z=1:size(train_coeff,2)
            temp1=squeeze(control_train_vep(x,y,:))';
            temp2=train_coeff(:,z);
            r=temp1*temp2;
            train_session_score(x,y,z)=r;
        end
    end
end

for x=1:size(control_test_vep,1)
    for y=1:size(control_test_vep,2)
        for z=1:size(train_coeff,2)
            temp1=squeeze(control_test_vep(x,y,:))';
            temp2=train_coeff(:,z);
            r=temp1*temp2;
            test_session_score(x,y,z)=r;
        end
    end
end

train_session1=squeeze(train_session_score(:,1,:));
train_session2=squeeze(train_session_score(:,2,:));

test_session1=squeeze(test_session_score(:,1,:));
test_session2=squeeze(test_session_score(:,2,:));

figure(4)
for x=1:PC_no
    subplot(ceil(PC_no/3),3,x)
    hold on
    plot(test_session1(:,x),test_session2(:,x),'.b','MarkerSize',12)
    plot([-2 2],[-2 2],'--k')
    train_corr=corrcoef(train_session1(:,x),train_session2(:,x));
    test_corr=corrcoef(test_session1(:,x),test_session2(:,x));
    ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-2.2 2.2];ax.YLim=[-2.2 2.2];
    axis('square')
    title(['PC' num2str(x) ' train corrcoef=' num2str(train_corr(1,2))])
    xlabel(['PC' num2str(x) ' test corrcoef=' num2str(test_corr(1,2))])
end


figure(5)
for x=1:size(test_session1,2)
    temp_corr=corrcoef(test_session1(:,x),test_session2(:,x));
    test_corr(x,:)=temp_corr(1,2);
    temp_corr=corrcoef(train_session1(:,x),train_session2(:,x));
    train_corr(x,:)=temp_corr(1,2);
end
hold on
plot(test_corr,'.b','MarkerSize',12)
plot([0 40],[0.7 0.7],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[1 40];ax.YLim=[-1 1];
axis('square')
xlabel('Principal component')
ylabel('Pearson correlation coefficient (session 1 vs. session 2)')


figure(7)
traintest_session1=cat(1,train_session1,test_session1);
traintest_session2=cat(1,train_session2,test_session2);
[Eu_dist]=calcEuclidean_dist(PC_no,test_session1,test_session2);
Eu_dist2=Eu_dist-diag(Eu_dist);
Eu_dist_test2(Eu_dist2==0)=NaN;
edges=-1:0.5:4; 
for x=1:size(train_score,1)
    subplot(8,5,x)
    histogram(Eu_dist2(x,:),edges,'Normalization','probability')
    inter_better_intra(x,:)=length(find(Eu_dist2(x,:)<0))./(length(Eu_dist2(x,:))-1);
    intra_rank(x,:)=length(find(Eu_dist2(x,:)<0))+1;
    hold on
    plot([0 0],[0 0.5],'--k')
     axis('square')
     ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-1 4];ax.YLim=[0 0.5];
    if x==38
        xlabel('Euclidean space for session 1 vs. session 2 across subjects relative to within subjects')
    end
    if x==15
        ylabel('Proportion of subjects')
    end
    if x==3
        title('Train and Test subjects')
    end
end

figure(8)
hold on
h=histogram(inter_better_intra.*100,[0:1:50],'Normalization','probability');
temp=h.Values;
plot([1 1],[0 0.8],'--k')
axis('square')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0 50];ax.YLim=[0 0.8];
ylabel('Proportion of subjects')
xlabel('Percent inter-subject sessions closer than intra-subject session')

figure(9)
edges=[1:1:10 80];
center=1.5:1:10.5;
hold on
h=histogram(intra_rank,edges,'Normalization','probability');
temp=h.Values;
for x=1:length(temp)
    sum_temp(x,:)=sum(temp(1:x));
end
plot(center,sum_temp,'-ok')
axis('square')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[1 11];ax.YLim=[0 1];
ylabel('Proportion of subjects')
xlabel('Rank order')


explained_test=var(test_score)./sum(var(test_score))*100;
disp(['variance explained by PCA model for test dataset:' num2str(sum(explained_test(1:PC_no)))]);

%% Determine peak amplitude and peak latency for test subjects

[peakanalysis_session1]=calcVEPpeak(xdata,squeeze(control_test_vep(:,1,:)));
[peakanalysis_session2]=calcVEPpeak(xdata,squeeze(control_test_vep(:,2,:)));



figure(15)

subplot(2,3,1)
hold on
test_PAcorr_N1L=corrcoef(peakanalysis_session1.N1_latency,peakanalysis_session2.N1_latency);
plot(peakanalysis_session1.N1_latency,peakanalysis_session2.N1_latency,'.b','MarkerSize',12)
plot([0.06 0.09],[0.06 0.09],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0.06 0.09];ax.YLim=[0.06 0.09];
axis('square')
xlabel('N1 latency (session 1)')
ylabel('N1 latency (session 2)')
title(['R=' num2str(test_PAcorr_N1L(1,2))])

subplot(2,3,2)
hold on
test_PAcorr_P1L=corrcoef(peakanalysis_session1.P1_latency,peakanalysis_session2.P1_latency);
plot(peakanalysis_session1.P1_latency,peakanalysis_session2.P1_latency,'.b','MarkerSize',12)
plot([0.09 0.13],[0.09 0.13],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0.09 0.13];ax.YLim=[0.09 0.13];
axis('square')
xlabel('P1 latency (session 1)')
ylabel('P1 latency (session 2)')
title(['R=' num2str(test_PAcorr_P1L(1,2))])


subplot(2,3,3)
hold on
test_PAcorr_N2L=corrcoef(peakanalysis_session1.N2_latency,peakanalysis_session2.N2_latency);
plot(peakanalysis_session1.N2_latency,peakanalysis_session2.N2_latency,'.b','MarkerSize',12)
plot([0.12 0.175],[0.12 0.175],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0.12 0.175];ax.YLim=[0.12 0.175];
axis('square')
xlabel('N2 latency (session 1)')
ylabel('N2 latency (session 2)')
title(['R=' num2str(test_PAcorr_N2L(1,2))])

subplot(2,3,4)
hold on
test_PAcorr_N1P1A=corrcoef(peakanalysis_session1.N1P1_amplitude,peakanalysis_session2.N1P1_amplitude);
plot(peakanalysis_session1.N1P1_amplitude,peakanalysis_session2.N1P1_amplitude,'.b','MarkerSize',12)
plot([0.005 0.5],[0.005 0.5],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0.005 0.5];ax.YLim=[0.005 0.5];
axis('square')
xlabel('N1P1 amplitude (session 1)')
ylabel('N1P1 amplitude (session 2)')
title(['R=' num2str(test_PAcorr_N1P1A(1,2))])

subplot(2,3,5)
hold on
test_PAcorr_P1N2A=corrcoef(peakanalysis_session1.P1N2_amplitude,peakanalysis_session2.P1N2_amplitude);
plot(peakanalysis_session1.P1N2_amplitude,peakanalysis_session2.P1N2_amplitude,'.b','MarkerSize',12)
plot([0.00003 0.3],[0.00003 0.3],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0.00003 0.3];ax.YLim=[0.00003 0.3];
axis('square')
xlabel('P1N2 amplitude (session 1)')
ylabel('P1N2 latency (session 2)')
title(['R=' num2str(test_PAcorr_P1N2A(1,2))])


%% Determine peak amplitude and peak latency for test subjects usng diopsys data

figure(16)

subplot(1,3,1)
hold on
test_PAcorr_N1L=corrcoef(control_vep_subjects_sess2.Left_Cursor_Lat,control_vep_subjects_sess22.Left_Cursor_Lat);
plot(control_vep_subjects_sess2.Left_Cursor_Lat,control_vep_subjects_sess22.Left_Cursor_Lat,'.b','MarkerSize',12)
plot([0.06 0.09],[0.06 0.09],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0 100];ax.YLim=[0 100];
axis('square')
xlabel('N1 latency (session 1)')
ylabel('N1 latency (session 2)')
title(['R=' num2str(test_PAcorr_N1L(1,2))])

subplot(1,3,2)
hold on
test_PAcorr_P1L=corrcoef(control_vep_subjects_sess2.Right_Cursor_Lat,control_vep_subjects_sess22.Right_Cursor_Lat);
plot(control_vep_subjects_sess2.Right_Cursor_Lat,control_vep_subjects_sess22.Right_Cursor_Lat,'.b','MarkerSize',12)
plot([0.09 0.13],[0.09 0.13],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0 130];ax.YLim=[0 130];
axis('square')
xlabel('P1 latency (session 1)')
ylabel('P1 latency (session 2)')
title(['R=' num2str(test_PAcorr_P1L(1,2))])


subplot(1,3,3)
hold on
test_PAcorr_N1P1A=corrcoef(control_vep_subjects_sess2.Delta_Amp,control_vep_subjects_sess22.Delta_Amp);
plot(control_vep_subjects_sess2.Delta_Amp,control_vep_subjects_sess22.Delta_Amp,'.b','MarkerSize',12)
plot([0.005 0.5],[0.005 0.5],'--k')
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[0 50];ax.YLim=[0 50];
axis('square')
xlabel('N1P1 amplitude (session 1)')
ylabel('N1P1 amplitude (session 2)')
title(['R=' num2str(test_PAcorr_N1P1A(1,2))])


