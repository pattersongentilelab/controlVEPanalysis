% Use single gamma model (can adjust to different number of gammas)


data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])
saveFig_path = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/';


%% split recording into two halves and compare metrics
control_subject = find(cleaned_vep_files.subjecttype=='Control');
cleaned_vep = cleaned_vep(control_subject,:);
control_vep_subjects = cleaned_vep_files(control_subject,:);
xdata = cleaned_vep{1,3}.*1000; %convert to ms

% Select first session for each participant
vep = zeros(1,2,length(xdata));
subject = control_vep_subjects(1,:);
uniqueID = unique(control_vep_subjects.uniqueID);

counter = 1;
for x = 1:length(uniqueID)
    temp_loc = find(cell2mat(cleaned_vep(:,1))==uniqueID(x));
    subject(x,:) = control_vep_subjects(temp_loc(1),:);
    temp = cleaned_vep{temp_loc(1),4};
    if size(temp,1)>150
        randAll = randperm(size(temp,1),size(temp,1)); % randomly split the data in two parts
        vep(counter,1,:) = mean(temp(randAll(1:round(end/2)),:),1).*100; %correct amplification issue with diopsys
        vep(counter,2,:) = mean(temp(randAll(round(end/2)+1:end),:),1).*100; %correct amplification issue with diopsys
        counter = counter+1;
    end
    clear temp_loc
end
clear uniqueID

% shift VEP to 0 = mean response
for x = 1:size(vep)
    for y = 1:2
        ydata = squeeze(vep(x,y,:));
        vep(x,y,:) = ydata - mean(ydata);
    end
end

%% Fit gamma model to selected subjects
% select time window and number of gammas in function
time_end = 500; % determines the epoch looked at in ms from t = 0 being the alternating checkerboard
xdata_end = length(xdata(xdata<=time_end));
gammaC = {'b','r','m','g','c','y'};
mdl_x = 0:0.1:time_end;
nGamma = 4;

xdata = xdata(:,1:xdata_end);
vep = vep(:,:,1:xdata_end);

meanVEP = squeeze(mean(vep,1));

r_val = zeros(size(vep,1),2);
r2 = zeros(size(vep,1),2);

%% Check for local minima and maxima to use in peak analysis

Amp75 = zeros(size(vep,1),2);
Peak75 = zeros(size(vep,1),2);
Amp100 = zeros(size(vep,1),2);
Peak100 = zeros(size(vep,1),2);
Amp135 = zeros(size(vep,1),2);
Peak135 = zeros(size(vep,1),2);
Amp220 = zeros(size(vep,1),2);
Peak220 = zeros(size(vep,1),2);

for i = 1:size(vep,1)
    for j = 1:2
        ydata = squeeze(vep(i,j,:));
        diffY = diff([min(ydata) max(ydata)]);
        min_loc = islocalmin(ydata,'MinProminence',diffY*.2);
        min_peak = xdata(min_loc==1);
        max_loc = islocalmax(ydata,'MinProminence',diffY*.2);
        max_peak = xdata(max_loc==1);
        
        x = sum(min_loc(xdata>60 & xdata<90));
        switch x
            case 0
                amp75 = min(ydata(xdata>60 & xdata<90));
                peak75 = xdata(ydata==amp75);
                peak75 = peak75(1);
            case 1
                peak75 = min_peak(min_peak>60 & min_peak<90);
                amp75 = ydata(xdata==peak75);
            otherwise
                peak75 = min_peak(min_peak>60 & min_peak<90);
                peak75 = peak75(1);
                amp75 = ydata(xdata==peak75);
        end
        if amp75>-1
            amp75 = -1;
        end
        
        
         x = sum(max_loc(xdata>90 & xdata<130));
        switch x
            case 0
                amp100 = max(ydata(xdata>90 & xdata<130));
                peak100 = xdata(ydata==amp100);
                peak100 = peak100(1);
            case 1
                peak100 = max_peak(max_peak>90 & max_peak<130);
                amp100 = ydata(xdata==peak100);
            otherwise
                peak100 = max_peak(max_peak>90 & max_peak<130);
                peak100 = peak100(1);
                amp100 = ydata(xdata==peak100);
        end
        
        if amp100<1
            amp100 = 1;
        end
       
        x = sum(min_loc(xdata>100 & xdata<200));
        switch x
            case 0
                amp135 = min(ydata(xdata>100 & xdata<200));
                peak135 = xdata(ydata==amp135);
                peak135 = peak135(1);
            case 1
                peak135 = min_peak(min_peak>100 & min_peak<200);
                amp135 = ydata(xdata==peak135);
            otherwise
                peak135 = min_peak(min_peak>100 & min_peak<200);
                peak135 = peak135(1);
                amp135 = ydata(xdata==peak135);
        end
        
        if amp135>-1
            amp135 = -1;
        end
        
       x = sum(max_loc(xdata>200 & xdata<350));
        switch x
            case 0
                amp220 = max(ydata(xdata>200 & xdata<350));
                peak220 = xdata(ydata==amp220);
                peak220 = peak220(1);
            case 1
                peak220 = max_peak(max_peak>200 & max_peak<350);
                amp220 = ydata(xdata==peak220(1));
            otherwise
                peak220 = max_peak(max_peak>200 & max_peak<350);
                peak220 = peak220(1);
                amp220 = ydata(xdata==peak220(1));
        end
 
        if amp220<1
            amp220 = 1;
        end
          
        amp(i,j,1) = amp75;
        peak(i,j,1) = peak75;
        amp(i,j,2) = amp100;
        peak(i,j,2) = peak100;
        amp(i,j,3) = amp135;
        peak(i,j,3) = peak135;
        amp(i,j,4) = amp220;
        peak(i,j,4) = peak220;
    end
end

%% Compare across split data

fig = figure(300);
z = 1;
for x = 1:2
    for y = 1:nGamma
        subplot(3,nGamma,z)
        hold on
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    
        switch x
            case 1
                plot(squeeze(peak(:,1,y)),squeeze(peak(:,2,y)),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})
                plot([0 400],[0 400],'--')
                [rr,pp] = corrcoef(squeeze(peak(:,1,y)),squeeze(peak(:,2,y)));
                title(sprintf('peak time %1d',y))
                ylim([0 400])
                xlim([0 400])
            case 2
                plot(squeeze(abs(amp(:,1,y))),squeeze(abs(amp(:,2,y))),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})
                plot([0 max(max(max(amp)))],[0 max(max(max(amp)))],'--')
                [rr,pp] = corrcoef(squeeze(abs(amp(:,1,y))),squeeze(abs(amp(:,2,y))));
                title(sprintf('amplitude %1d',y))
                ylim([0 max(max(max(amp)))])
                xlim([0 max(max(max(amp)))])
        end
        xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
        z = z+1;
    end
end


%% Bland-Altman plot

% Informed basis set

fig = figure(302);

subplot(2,4,1)
diff_n75pk_sess = diff(squeeze(peak(:,:,1)),[],2);
mean_n75pk_sess = mean(squeeze(peak(:,:,1)),2);
hold on
plot(mean_n75pk_sess,diff_n75pk_sess,'o','MarkerFaceColor','b','MarkerEdgeColor','w','Color','b')
plot([min(mean_n75pk_sess) max(mean_n75pk_sess)],[mean(diff_n75pk_sess) mean(diff_n75pk_sess)],'--k')
plot([min(mean_n75pk_sess) max(mean_n75pk_sess)],[mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96 mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96],'--k')
plot([min(mean_n75pk_sess) max(mean_n75pk_sess)],[mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96 mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96],'--k')
text(max(mean_n75pk_sess)-2,mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96));
text(max(mean_n75pk_sess)-2,mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96+2,sprintf('%2.1f',mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96));
title('Bland Altman N75 peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75pk_sess), max(mean_n75pk_sess)]; ax.YLim = [min(diff_n75pk_sess)-1,max(diff_n75pk_sess)+1];

subplot(2,4,2)
diff_p100pk_sess = diff(squeeze(peak(:,:,2)),[],2);
mean_p100pk_sess = mean(squeeze(peak(:,:,2)),2);
hold on
plot(mean_p100pk_sess,diff_p100pk_sess,'o','MarkerFaceColor','r','MarkerEdgeColor','w','Color','r')
plot([min(mean_p100pk_sess) max(mean_p100pk_sess)],[mean(diff_p100pk_sess) mean(diff_p100pk_sess)],'--k')
plot([min(mean_p100pk_sess) max(mean_p100pk_sess)],[mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96 mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96],'--k')
plot([min(mean_p100pk_sess) max(mean_p100pk_sess)],[mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96 mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96],'--k')
text(max(mean_p100pk_sess)-2,mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96));
text(max(mean_p100pk_sess)-2,mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96+2,sprintf('%2.1f',mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96));
title('Bland Altman P100 peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_p100pk_sess), max(mean_p100pk_sess)]; ax.YLim = [min(diff_p100pk_sess)-1,max(diff_p100pk_sess)+1];

subplot(2,4,3)
diff_n135pk_sess = diff(squeeze(peak(:,:,3)),[],2);
mean_n135pk_sess = mean(squeeze(peak(:,:,3)),2);
hold on
plot(mean_n135pk_sess,diff_n135pk_sess,'o','MarkerFaceColor','m','MarkerEdgeColor','w','Color','m')
plot([min(mean_n135pk_sess) max(mean_n135pk_sess)],[mean(diff_n135pk_sess) mean(diff_n135pk_sess)],'--k')
plot([min(mean_n135pk_sess) max(mean_n135pk_sess)],[mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96 mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96],'--k')
plot([min(mean_n135pk_sess) max(mean_n135pk_sess)],[mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96 mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96],'--k')
text(max(mean_n135pk_sess)-2,mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96));
text(max(mean_n135pk_sess)-2,mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96+2,sprintf('%2.1f',mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96));
title('Bland Altman N135 peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n135pk_sess), max(mean_n135pk_sess)]; ax.YLim = [min(diff_n135pk_sess)-1,max(diff_n135pk_sess)+1];

subplot(2,4,4)
diff_pLatepk_sess = diff(squeeze(peak(:,:,4)),[],2);
mean_pLatepk_sess = mean(squeeze(peak(:,:,4)),2);
hold on
plot(mean_pLatepk_sess,diff_pLatepk_sess,'o','MarkerFaceColor','g','MarkerEdgeColor','w','Color','g')
plot([min(mean_pLatepk_sess) max(mean_pLatepk_sess)],[mean(diff_pLatepk_sess) mean(diff_pLatepk_sess)],'--k')
plot([min(mean_pLatepk_sess) max(mean_pLatepk_sess)],[mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96 mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96],'--k')
plot([min(mean_pLatepk_sess) max(mean_pLatepk_sess)],[mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96 mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96],'--k')
text(max(mean_pLatepk_sess)-2,mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96));
text(max(mean_pLatepk_sess)-2,mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96+2,sprintf('%2.1f',mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96));
title('Bland Altman Late peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_pLatepk_sess), max(mean_pLatepk_sess)]; ax.YLim = [min(diff_pLatepk_sess)-1,max(diff_pLatepk_sess)+1];


subplot(2,4,5)
diff_n75am_sess = diff(squeeze(amp(:,:,1)),[],2);
mean_n75am_sess = mean(squeeze(amp(:,:,1)),2);
hold on
plot(mean_n75am_sess,diff_n75am_sess,'o','MarkerFaceColor','b','MarkerEdgeColor','w','Color','b')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess) mean(diff_n75am_sess)],'--k')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96 mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96],'--k')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96 mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96],'--k')
text(max(mean_n75am_sess)-2,mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96));
text(max(mean_n75am_sess)-2,mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96+2,sprintf('%2.1f',mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96));
title('Bland Altman N75 amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75am_sess), max(mean_n75am_sess)]; ax.YLim = [min(diff_n75am_sess)-1,max(diff_n75am_sess)+1];

subplot(2,4,6)
diff_p100am_sess = diff(squeeze(amp(:,:,2)),[],2);
mean_p100am_sess = mean(squeeze(amp(:,:,2)),2);
hold on
plot(mean_p100am_sess,diff_p100am_sess,'o','MarkerFaceColor','r','MarkerEdgeColor','w','Color','r')
plot([min(mean_p100am_sess) max(mean_p100am_sess)],[mean(diff_p100am_sess) mean(diff_p100am_sess)],'--k')
plot([min(mean_p100am_sess) max(mean_p100am_sess)],[mean(diff_p100am_sess)+std(diff_p100am_sess)*1.96 mean(diff_p100am_sess)+std(diff_p100am_sess)*1.96],'--k')
plot([min(mean_p100am_sess) max(mean_p100am_sess)],[mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96 mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96],'--k')
text(max(mean_p100am_sess)-2,mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96-2,sprintf('%2.1f',mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96));
text(max(mean_p100am_sess)-2,mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96+2,sprintf('%2.1f',mean(diff_p100am_sess)+std(diff_p100am_sess)*1.96));
title('Bland Altman P100 amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_p100am_sess), max(mean_p100am_sess)]; ax.YLim = [min(diff_p100am_sess)-1,max(diff_p100am_sess)+1];

subplot(2,4,7)
diff_n135am_sess = diff(squeeze(amp(:,:,3)),[],2);
mean_n135am_sess = mean(squeeze(amp(:,:,3)),2);
hold on
plot(mean_n135am_sess,diff_n135am_sess,'o','MarkerFaceColor','m','MarkerEdgeColor','w','Color','m')
plot([min(mean_n135am_sess) max(mean_n135am_sess)],[mean(diff_n135am_sess) mean(diff_n135am_sess)],'--k')
plot([min(mean_n135am_sess) max(mean_n135am_sess)],[mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96 mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96],'--k')
plot([min(mean_n135am_sess) max(mean_n135am_sess)],[mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96 mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96],'--k')
text(max(mean_n135am_sess)-2,mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96));
text(max(mean_n135am_sess)-2,mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96+2,sprintf('%2.1f',mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96));
title('Bland Altman N135 amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n135am_sess), max(mean_n135am_sess)]; ax.YLim = [min(diff_n135am_sess)-1,max(diff_n135am_sess)+1];

subplot(2,4,8)
diff_pLateam_sess = diff(squeeze(amp(:,:,4)),[],2);
mean_pLateam_sess = mean(squeeze(amp(:,:,4)),2);
hold on
plot(mean_pLateam_sess,diff_pLateam_sess,'o','MarkerFaceColor','g','MarkerEdgeColor','w','Color','g')
plot([min(mean_pLateam_sess) max(mean_pLateam_sess)],[mean(diff_pLateam_sess) mean(diff_pLateam_sess)],'--k')
plot([min(mean_pLateam_sess) max(mean_pLateam_sess)],[mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96 mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96],'--k')
plot([min(mean_pLateam_sess) max(mean_pLateam_sess)],[mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96 mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96],'--k')
text(max(mean_pLateam_sess)-2,mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96-2,sprintf('%2.1f',mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96));
text(max(mean_pLateam_sess)-2,mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96+2,sprintf('%2.1f',mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96));
title('Bland Altman Late amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_pLateam_sess), max(mean_pLateam_sess)]; ax.YLim = [min(diff_pLateam_sess)-1,max(diff_pLateam_sess)+1];
