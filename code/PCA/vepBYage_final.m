% Fit control subject data to PCA model

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
analysis_path = getpref('controlVEPanalysis','controlVEP_AnalysisPath');
load([data_path '/cleaned_VEP.mat'])

xdata=cell2mat(cleaned_vep(1,3));
PC_no=7;
edges=-2.5:0.25:2.5;
center=edges(1:end-1)+diff(edges);
X=-2.5:0.01:2.5;
neuro_active = 0; % 0 = includes subjects on neuro-active meds, 1 = excludes subjects on
% neuroactive meds


%% Subject selection
control_subject=find(cleaned_vep_files.subjecttype=='Control');
cleaned_vep=cleaned_vep(control_subject,:);
control_vep_subject_sessions=cleaned_vep_files(control_subject,:);

if neuro_active==1
    load([data_path '/neuroActive_meds.mat']) 
    counter=1;
    for x=1:length(control_subject)
        Control_subject(counter,:)=control_subject(x);
        if isempty(find(neuro_active_meds==control_vep_subject_sessions.uniqueID(x)))==1
            control_subject2(counter,:)=x;
            counter=counter+1;
        end
    end
    cleaned_vep=cleaned_vep(control_subject2,:);
    control_vep_subject_sessions=control_vep_subject_sessions(control_subject2,:);
    control_subject=control_subject(control_subject2);
end

unique_ID=unique(control_vep_subject_sessions.uniqueID);

%% Load PCA model
load([analysis_path '/PCAmodel.mat'])


%% Fit PCA model to the whole dataset

TEMP=[];

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
        TEMP2=cell2mat(cleaned_vep(temp_loc(y),4));
        temp_ydata2(y,:)=nanmean(cell2mat(temp_ydata(y,:)));
        TEMP=cat(1,TEMP,TEMP2);
    end
    control_vep{x,:}=TEMP;
    TEMP=[];
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

age=control_vep_subjects.age_vep;


%% ANOVA for Age
mdl_Age=fitlm(controlM_score(:,1:7),age);

anova(mdl_Age,'summary')

% fit PC2 and PC3
X=controlM_score(:,2);
Y=controlM_score(:,3);

XY=controlM_score(:,[2 4]);
r=mean(XY);
[~,~,V]=svd(XY,0);
x_fit=@(z_fit) r(1)+(z_fit-r(3))/V(3,1)*V(1,1);
y_fit=@(z_fit) r(2)+(z_fit-r(3))/V(3,1)*V(2,1);


%% 2D scatter plot of PC2 and PC3 by age


% Open a figure
fig5 = figure(5);
fig5.Renderer='Painters';
subplot(5,4,1:16);

% What alpha do we want for the scatter?
scatterAlpha = 0.5;

% define color by age
minAge=10;
maxAge=20;
crange=[minAge maxAge];

cmap1=cat(1,fliplr(0:0.5/128:0.5),fliplr(0:0.5/128:0.5),0.5:0.5/128:1)';
cmap2=cat(1,0:0.5/128:0.5,0:0.5/128:0.5,0:0.5/128:0.5)';
cmap=cat(1,cmap2(2:end,:),cmap1(2:end,:));

% Normalize the values to be between 1 and 256
age(age < crange(1)) = crange(1);
age(age > crange(2)) = crange(2);
ageN = round(((age - crange(1)) ./ diff(crange)) .* 255)+1;
% Convert any nans to ones
ageN(isnan(ageN)) = 1;
% Convert the normalized values to the RGB values of the colormap
rgb = cmap(ageN, :);    
for ii = 1:length(age)
        % Normalize the values to be between 1 and 256 for cell ii
        ageN = age(ii);
        ageN(ageN < crange(1)) = crange(1);
        ageN(ageN > crange(2)) = crange(2);
        ageN = round(((ageN - crange(1)) ./ diff(crange)) .* 255)+1;
        % Convert any nans to ones
        ageN(isnan(ageN)) = 1;
        % Convert the normalized values to the RGB values of the colormap
        rgb(ii,:) = cmap(ageN, :);
end
    
C=rgb;

% Render the scatter plot; set the alpha
sHandle = scatter(X,Y,200,C,'filled');
sHandle.MarkerFaceAlpha = scatterAlpha;
sHandle.MarkerEdgeColor = [0.5 0.5 0.5];

% Adjust the limits, add labels, set the axes equal
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-2 2];ax.YLim=[-1 1];
xlabel('PC2')
ylabel('PC3')
axis equal

% We flip the direction of the X axis so that age increases left to right
set ( gca, 'xdir', 'reverse' )

% Find the coefficient values for the vector that reflects the effect of
% age
simAges = [10,20];
pcDims = [2, 3, 6];
xyz = nan(3);
uvw = nan(3);
for ii = 1:length(pcDims)
    temp=polyfit(age,controlM_score(:,pcDims(ii)),1);
    xyz(ii) = polyval(temp,simAges(1));
    uvw(ii) = polyval(temp,simAges(2))-polyval(temp,simAges(1));
end
hold on

% Draw an arrow for this vector
plot3([xyz(1) uvw(1)],[xyz(2) uvw(2)],[xyz(3) uvw(3)],'-k','LineWidth',2);

% Create a colorbar
subplot(5,4,17:20);
im = repmat(linspace(0,254,100),10,1);
colormap(cmap);
h = image(im);
h.AlphaData = scatterAlpha;
nTicks = 6;
xticks(linspace(1,100,nTicks));
xticklabels(linspace(crange(1),crange(2),nTicks));
h = gca; h.YAxis.Visible = 'off';
box off

%% compare age
figure(10)
for x=1:PC_no
    subplot(ceil(PC_no/3),3,x)
    hold on
    plot(control_vep_subjects.age_vep,controlM_score(:,x),'.','Color',[0.5 0.5 0.5],'MarkerSize',12)
    ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[-2 2];
    axis('square')
    corr=corrcoef(control_vep_subjects.age_vep,controlM_score(:,x));
    lsline
    title(['PC' num2str(x) ' corrcoef=' num2str(corr(1,2))])
    if x==PC_no 
        xlabel('Age at time of first VEP')
    end
    if x==4
        ylabel('Component score')
    end
end


%% Age effects on peak analysis

[peak_analysis]=calcVEPpeak(xdata,controlM_vep);

figure(30)

subplot(2,3,1)
hold on
plot(mAge,peak_analysis.N1_latency,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(peak_analysis.N1_latency) max(peak_analysis.N1_latency)];
axis('square')
ylabel('N1 latency')
xlabel('age')
corr_age=corrcoef(mAge,peak_analysis.N1_latency);
title(['R=' num2str(corr_age(1,2))])

subplot(2,3,2)
hold on
plot(mAge,peak_analysis.P1_latency,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(peak_analysis.P1_latency) max(peak_analysis.P1_latency)];
axis('square')
ylabel('P1 latency')
xlabel('age')
corr_age=corrcoef(mAge,peak_analysis.P1_latency);
title(['R=' num2str(corr_age(1,2))])

subplot(2,3,3)
hold on
plot(mAge,peak_analysis.N2_latency,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(peak_analysis.N2_latency) max(peak_analysis.N2_latency)];
axis('square')
ylabel('N2 latency')
xlabel('age')
corr_age=corrcoef(mAge,peak_analysis.N2_latency);
title(['R=' num2str(corr_age(1,2))])

subplot(2,3,4)
hold on
plot(mAge,peak_analysis.N1P1_amplitude,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(peak_analysis.N1P1_amplitude) max(peak_analysis.N1P1_amplitude)];
axis('square')
xlabel('age')
ylabel('N1P1 amplitude')
corr_age=corrcoef(mAge,peak_analysis.N1P1_amplitude);
title(['R=' num2str(corr_age(1,2))])

subplot(2,3,5)
hold on
plot(mAge,peak_analysis.P1N2_amplitude,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(peak_analysis.P1N2_amplitude) max(peak_analysis.P1N2_amplitude)];
axis('square')
xlabel('age')
ylabel('P2N1 amplitude')
corr_age=corrcoef(mAge,peak_analysis.P1N2_amplitude);
title(['R=' num2str(corr_age(1,2))])


%% Age effects on peak analysis using diopsys data


figure(31)

subplot(1,3,1)
hold on
plot(mAge,control_vep_subjects.Left_Cursor_Lat,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(control_vep_subjects.Left_Cursor_Lat) max(control_vep_subjects.Left_Cursor_Lat)];
axis('square')
ylabel('N1 latency')
xlabel('age')
corr_age=corrcoef(mAge,control_vep_subjects.Left_Cursor_Lat);
title(['R=' num2str(corr_age(1,2))])

subplot(1,3,2)
hold on
plot(mAge,control_vep_subjects.Right_Cursor_Lat,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(control_vep_subjects.Right_Cursor_Lat) max(control_vep_subjects.Right_Cursor_Lat)];
axis('square')
ylabel('P1 latency')
xlabel('age')
corr_age=corrcoef(mAge,control_vep_subjects.Right_Cursor_Lat);
title(['R=' num2str(corr_age(1,2))])

subplot(1,3,3)
hold on
plot(mAge,control_vep_subjects.Delta_Amp,'.b','MarkerSize',12)
lsline
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[min(control_vep_subjects.Delta_Amp) max(control_vep_subjects.Delta_Amp)];
axis('square')
xlabel('age')
ylabel('N1P1 amplitude')
corr_age=corrcoef(mAge,control_vep_subjects.Delta_Amp);
title(['R=' num2str(corr_age(1,2))])


%% model based on head circumference
hc=control_vep_subjects.head_circumference;
hc(hc<0)=NaN;

has_hc=find(isnan(hc)==0);

figure(30)
hold on
plot(age(has_hc),hc(has_hc),'ob')
corr=corrcoef(age(has_hc),hc(has_hc));
ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[10 20];ax.YLim=[18 25];
axis('square')
title(['R=' num2str(corr(1,2))])
xlabel('Age at VEP')
ylabel('Head circumference')
    
mdl_HC=fitlm(controlM_score(:,1:7),hc);

anova(mdl_HC,'summary')

%% regress out effect of HC and retest age model

[ADJcontrolM_score]=AgeSexHcNull(PC_no,controlM_score(has_hc,1:7),[],[],hc(has_hc));

mdl_Age_hc=fitlm(ADJcontrolM_score(:,1:7),age(has_hc));

anova(mdl_Age_hc,'summary')

%% Simulated VEPs

% Simulate a 10 year old, 15 year old, and 20 year old sex-nulled subject
[simVEP]=synthVEP(7,controlM_score,train_coeff,explained','age',age);


%% Sex comparison
sex=control_vep_subjects.sex_master;

mdl_Sex=fitlm(controlM_score(:,1:7),sex);

disp('sex unadjusted')
anova(mdl_Sex,'summary')


[ADJAcontrolM_score]=AgeSexNull(PC_no,controlM_score,[],age);

mdl_Sex2=fitlm(ADJAcontrolM_score,sex);

disp('sex, age adjusted')
anova(mdl_Sex2,'summary')

female_subjects=find(control_vep_subjects.sex_master==1);
male_subjects=find(control_vep_subjects.sex_master==2);


 [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,controlM_vep(female_subjects,:),...
                'errorbars','Boot','color_mean',[0 0 0],'color_err',[0.5 0.5 0.5],'fig_num',61,...
                'sub_plot',true,'sub_plot_num',[1 3 1]);
        
  [y_dataM,y_dataERR1,y_dataERR2]=plot_meanVEP(xdata,controlM_vep(male_subjects,:),...
                'errorbars','Boot','color_mean',[0 1 1],'color_err',[0.5 1 1],'fig_num',61,...
                'sub_plot',true,'sub_plot_num',[1 3 1]);
            
  legend('show')
  legend('females','','males','')
  
  

%% compare sexes
female_subjects=find(control_vep_subjects.sex_master==1);
male_subjects=find(control_vep_subjects.sex_master==2);

figure(11)
edges=-2.5:0.25:2.5;
center=edges(1:end-1)+diff(edges);
X=-2.5:0.01:2.5;
for x=1:PC_no
    subplot(ceil(PC_no/3),3,x)
    hold on
    hist1=histcounts(controlM_score(female_subjects,x),edges,'Normalization','probability');
    Y1=spline(center,hist1,X);
    hist2=histcounts(controlM_score(male_subjects,x),edges,'Normalization','probability');
    Y2=spline(center,hist2,X);
    plot(center,hist1,'.k',X,Y1,'-k')
    plot(center,hist2,'.b',X,Y2,'-b')
    [h p]=kstest2(controlM_score(female_subjects,x),controlM_score(male_subjects,x));
    ax=gca;ax.Box='off';ax.TickDir='out';ax.XLim=[-2.5 2.5];ax.YLim=[0 0.5];
    title(['p-value=' num2str(p)]);
end
