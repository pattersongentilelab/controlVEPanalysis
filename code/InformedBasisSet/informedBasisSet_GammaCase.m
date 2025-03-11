% Use single gamma model (can adjust to different number of gammas)

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])



case_subject = find(cleaned_vep_files.subjecttype=='Case');
cleaned_vep = cleaned_vep(case_subject,:);
case_vep_subjects = cleaned_vep_files(case_subject,:);
xdata = cleaned_vep{1,3}.*1000; %convert to ms

case_vep_subjects.migrainePCSI = sum([case_vep_subjects.teenPCSI_headache case_vep_subjects.teenPCSI_senslight case_vep_subjects.teenPCSI_sensnoise case_vep_subjects.teenPCSI_nausea],2);

unique_ID = unique(case_vep_subjects.uniqueID);

%% select subjects who were seen within 28 days of concussion

% Select first session for each participant
vep = zeros(1,length(xdata));
subject = case_vep_subjects(1,:);
uniqueID = unique(case_vep_subjects.uniqueID);

signal_var = NaN*ones(69,1);
noise_var = NaN*ones(69,1);
counter = 1;
for x = 1:length(uniqueID)
    temp_loc = find(cell2mat(cleaned_vep(:,1))==uniqueID(x));
    if size(cleaned_vep{temp_loc(1),4},1)>150
        vep(counter,:) = mean(cleaned_vep{temp_loc(1),4},1);
        subject(counter,:) = case_vep_subjects(temp_loc(1),:);
        temp = cleaned_vep{temp_loc(1),4};
        temp2 = mean(cleaned_vep{temp_loc(1),4},1);
        signal_var(counter,:) = var(temp2);
        noise_var(counter,:) = mean(var(temp-temp2));
        counter = counter+1;
    end
    clear temp_loc
end
clear uniqueID

subject.snr = signal_var./noise_var;

% filled out teen PCSI and was seen within 28 days of concussion
subject = subject(subject.dayspostinj<=28 & ~isnan(subject.teenPCSI_score_total),:);
vep = vep(subject.dayspostinj<=28 & ~isnan(subject.teenPCSI_score_total),:);

% determine days post injury symptom free

subject.conc_date = subject.VEPraw_date(:,1) - days(subject.dayspostinj(:,1));

% for subjects without a recovery date, look at the last PCSI - if it is
% <7, that will be counted as recovery date

for i = 1:height(subject)
    if isnat(subject.date_sympfree(i))
        temp = case_vep_subjects(case_vep_subjects.uniqueID==subject.uniqueID(i),:);
        temp2 = temp(temp.teenPCSI_score_total<7,:);
        if ~isempty(temp2)
            subject.date_sympfree(i) = temp2.pcsidate(1);
        end
    end
end

subject.dayspost_recov = between(subject.conc_date,subject.date_sympfree,'Days');
subject.dayspost_recov = split(subject.dayspost_recov,'d');

for i = 1:height(subject)
    if isnan(subject.dayspost_recov(i))
        temp = case_vep_subjects(case_vep_subjects.uniqueID==subject.uniqueID(i),:);
        if temp.dayspostinj>28
            subject.dayspost_recov(i) = -1;
        end
    end
end

subject.persist28 = NaN*ones(height(subject),1);
subject.persist28(subject.dayspost_recov<=28) = 0;
subject.persist28(subject.dayspost_recov>28) = 1;

% look for PCSI score >28 days, or if they recovered in < 28 days use the last PCSI
% score recorded

for j = 1:height(subject)
    temp = case_vep_subjects(case_vep_subjects.uniqueID==subject.uniqueID(j) & ~isnan(case_vep_subjects.teenPCSI_score_total),:);
    switch subject.persist28(j)
        case 0
            subjectFU(j,:) = temp(end,:);
        case 1
            temp2 = temp(temp.dayspostinj>28,:);
            if ~isempty(temp2)  
                subjectFU(j,:) = temp2(1,:);
            end
    end
end

 subject = subject(subjectFU.uniqueID~=0,:);
 subjectFU = subjectFU(subjectFU.uniqueID~=0,:);
 vep = vep(subjectFU.uniqueID~=0,:); 
 
 subjectFU.dayspost_recov = subject.dayspost_recov;
 subjectFU.persist28 = subject.persist28;
 subjectFU.snr =  subject.snr;
 
subject.Mig = zeros(height(subject),1);
subject.Mig(((subject.teenPCSI_senslight>0 & subject.teenPCSI_sensnoise>0)|subject.teenPCSI_nausea>0) & subject.teenPCSI_headache>=3) = 1;
subjectFU.Mig = zeros(height(subjectFU),1);
subjectFU.Mig(((subjectFU.teenPCSI_senslight>0 & subjectFU.teenPCSI_sensnoise>0)|subjectFU.teenPCSI_nausea>0) & subjectFU.teenPCSI_headache>=3) = 1;

% shift VEP to 0 = mean response
for x = 1:size(vep)
        ydata = vep(x,:);
        ydata = ydata - mean(ydata);
        ydata = ydata./max(abs(ydata)); % normalize by absolute maximum response
        vep(x,:) = ydata;
end

%% Fit gamma model to selected subjects
% select time window and number of gammas in function
time_end = 500; % determines the epoch looked at in ms from t = 0 being the alternating checkerboard
xdata_end = length(xdata(xdata<=time_end));
gammaC = {'b','r','m','g','c','y'};
mdl_x = 0:0.1:time_end;
nGamma = 4;

xdata = xdata(1:xdata_end);
vep = vep(:,1:xdata_end);

meanVEP = mean(vep,1);

r_val = zeros(size(vep,1),1);
r_val300 = zeros(size(vep,1),1);


%% Determine fits on individual VEP data
mdl = zeros(size(vep,1),3*nGamma);
gamma = zeros(nGamma,length(mdl_x));
Gamma = zeros(size(vep,1),nGamma,length(mdl_x));
yFit = zeros(size(vep));
bandwidth = NaN*ones(size(vep,1),nGamma);

Amp75 = zeros(size(vep,1),1);
Peak75 = zeros(size(vep,1),1);
Bw75 = zeros(size(vep,1),1);
Amp100 = zeros(size(vep,1),1);
Peak100 = zeros(size(vep,1),1);
Bw100 = zeros(size(vep,1),1);
Amp135 = zeros(size(vep,1),1);
Peak135 = zeros(size(vep,1),1);
Bw135 = zeros(size(vep,1),1);
Amp220 = zeros(size(vep,1),1);
Peak220 = zeros(size(vep,1),1);
Bw220 = zeros(size(vep,1),1);

for i = 1:size(vep,1)
        ydata = vep(i,:);
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
        if amp75>-0.01
            amp75 = -0.01;
        end
                
         x = sum(max_loc(xdata>peak75+5 & xdata<130));
        switch x
            case 0
                amp100 = max(ydata(xdata>peak75+5 & xdata<130));
                peak100 = xdata(ydata==amp100);
                peak100 = peak100(1);
            case 1
                peak100 = max_peak(max_peak>peak75+5 & max_peak<130);
                amp100 = ydata(xdata==peak100);
            otherwise
                peak100 = max_peak(max_peak>peak75+5 & max_peak<130);
                peak100 = peak100(1);
                amp100 = ydata(xdata==peak100);
        end
        
        if amp100<0.01
            amp100 = 0.01;
        end
       
        x = sum(min_loc(xdata>peak100+5 & xdata<200));
        switch x
            case 0
                amp135 = min(ydata(xdata>peak100+5 & xdata<200));
                peak135 = xdata(ydata==amp135);
                peak135 = peak135(1);
            case 1
                peak135 = min_peak(min_peak>peak100+5 & min_peak<200);
                amp135 = ydata(xdata==peak135);
            otherwise
                peak135 = min_peak(min_peak>peak100+5 & min_peak<200);
                peak135 = peak135(1);
                amp135 = ydata(xdata==peak135);
        end
        
        if amp135>-0.01
            amp135 = -0.01;
        end
        
       x = sum(max_loc(xdata>peak135+30 & xdata<350));
        switch x
            case 0
                amp220 = max(ydata(xdata>peak135+30 & xdata<350));
                peak220 = xdata(ydata==amp220);
                peak220 = peak220(1);
            case 1
                peak220 = max_peak(max_peak>peak135+30 & max_peak<350);
                amp220 = ydata(xdata==peak220);
            otherwise
                peak220 = max_peak(max_peak>peak135+30 & max_peak<350);
                peak220 = peak220(1);
                amp220 = ydata(xdata==peak220);
        end
 
        if amp220<0.01
            amp220 = 0.01;
        end
        
        bw75 = 10^((80-abs(diff([peak75 peak100])))/30);

        if bw75 < 30
            bw75 = 30;
        end
        if bw75 > 80
            bw75 = 80;
        end
        bw100 = 10^((100-(0.5*abs(diff([peak75 peak135]))))/40);
        if bw100 < 20
            bw100 = 20;
        end
        if bw100 > 80
            bw100 = 80;
        end
        bw135 = 10^((150-(0.5*abs(diff([peak100 peak220]))))/60);
        if bw135 < 20
            bw135 = 20;
        end
        if bw135 > 80
            bw135 = 80;
        end
        bw220 = 10^((200-abs(diff([peak135 peak220])))/80);
        if bw220 < 10
            bw220 = 10;
        end
        if bw220 > 80
            bw220 = 80;
        end
                
        Amp75(i,:) = amp75;
        Peak75(i,:) = peak75;
        Bw75(i,:) = bw75;
        Amp100(i,:) = amp100;
        Peak100(i,:) = peak100;
        Bw100(i,:) = bw100;
        Amp135(i,:) = amp135;
        Peak135(i,:) = peak135;
        Bw135(i,:) = bw135;
        Amp220(i,:) = amp220;
        Peak220(i,:) = peak220;
        Bw220(i,:) = bw220;

end

for i = 1:size(vep,1)

    ydata = vep(i,:);
    
    p0 = [Bw75(i,:) Peak75(i,:) Amp75(i,:) Bw100(i,:) Peak100(i,:) Amp100(i,:) Bw135(i,:) Peak135(i,:) Amp135(i,:) Bw220(i,:) Peak220(i,:) Amp220(i,:)];
    lb = [max([30 Bw75(i,:)-5]) Peak75(i,:)-2 Amp75(i,:)*1.1 max([20 Bw100(i,:)-5]) Peak100(i,:)-3 Amp100(i,:)*0.9 max([15 Bw135(i,:)-5]) Peak135(i,:)-5 Amp135(i,:)*1.1 max([15 Bw220(i,:)-5]) Peak220(i,:)-5 Amp220(i,:)*0.9]; 
    ub = [min([110 Bw75(i,:)+5]) Peak75(i,:)+2 Amp75(i,:)*0.9 min([110 Bw100(i,:)+5]) Peak100(i,:)+3 Amp100(i,:)*1.1 min([110 Bw135(i,:)+5]) Peak135(i,:)+5 -Amp135(i,:)*0.9 min([100 Bw220(i,:)+5]) Peak220(i,:)+5 Amp220(i,:)*1.1];

    myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,nGamma)).^2));
    mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [vep_fit,gamma] = gammaVEP_model(mdl_x,mdl(i,:),nGamma);
    Gamma(i,:,:) = gamma;
    
    for z = 1:nGamma
        bandwidth(i,z,:) = gamma_bandwidth(mdl_x,gamma(z,:));
    end

    [yFit(i,:)] = gammaVEP_model(xdata,mdl(i,:),nGamma);
    r = corrcoef(ydata,yFit(i,:));
    r_val(i,:) = r(1,2);
    r = corrcoef(ydata(1,1:307),yFit(i,1:307));
    r_val300(i,:) = r(1,2);

    switch i
        case {1,9,17,25,33,41,49,57}
            figure
            j = 1;
        case {5,13,21,29,37,45,53}
            j = 9;
        otherwise
            j = j+1;
    end

    subplot(4,4,j)
    plot(xdata,vep(i,:),'-','Color',[0.5 0.5 0.5])
    hold on
    plot(xdata,yFit(i,:),'-k','LineWidth',2)
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 time_end]; ax.YLim = [-2 2];
    xlabel(sprintf('r = %2.2f',r_val(i)))
    title(sprintf('%2.0f, C%2.0f, age = %2.0f',[i subject.uniqueID(i) subject.age_vep(i)]))
    subplot(4,4,j+4)
    hold on
    for X = 1:nGamma
         plot(mdl_x,gamma(X,:),['-' gammaC{X}])
    end
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-2 2]; ax.XLim = [0 time_end];
    
    switch i
        case {8,16,24,32,40,48,56}
    end
end

peak = mdl(:,2:3:end);
amp = mdl(:,3:3:end);

    
% Plot R values
figure(202)
hold on
plot(ones(size(r_val)),r_val,'.','Color',[0.5 0.5 0.5])
errorbar(1,mean(r_val),std(r_val),'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
plot(2*ones(size(r_val300)),r_val300,'.','Color',[0.5 0.5 0.5])
errorbar(2,mean(r_val300),std(r_val300),'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 3];



%% Plot mean VEP and parameters


figure(205)
subplot(1,2,1)
hold on
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';

figure(205)
subplot(1,2,2)
hold on
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 5]; title('recovered within 28 days');

% plotting parameters

for y = 1:nGamma
     
    plot_meanVEP(mdl_x,squeeze(Gamma(subject.Mig==0,y,:)),'errorbars','Boot','fig_num',205,'sub_plot',true,'sub_plot_num',[1 2 1],'color_mean',[0 0 0],'color_err',[0.8 0.8 0.8]);
    
    figure(205)
    subplot(1,2,1)
    bootPk = bootstrp(1000,@mean,peak(subject.Mig==0,y));
    bootPk = sort(bootPk,1);
    bootAmp = bootstrp(1000,@mean,amp(subject.Mig==0,y));
    bootAmp = sort(bootAmp,1);
    errorbar(bootPk(500),bootAmp(500),abs(diff(bootAmp([25 500]))),abs(diff(bootAmp([500 975]))),abs(diff(bootPk([25 500]))),abs(diff(bootPk([500 975]))),'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','w','Color','k')
    
    figure(205)
    subplot(1,2,2)
    hold on
    bootBW = bootstrp(1000,@mean,bandwidth(subject.Mig==0,y));
    bootBW = sort(bootBW,1);
    errorbar(y,bootBW(500),abs(diff(bootBW([25 500]))),abs(diff(bootBW([500 975]))),'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','w','Color','k')

end




for y = 1:nGamma
    plot_meanVEP(mdl_x,squeeze(Gamma(subject.Mig==1,y,:)),'errorbars','Boot','fig_num',205,'sub_plot',true,'sub_plot_num',[1 2 1],'color_mean',[1 0 0],'color_err',[1 0.8 0.8]);
    figure(205)
    subplot(1,2,1)
    bootPk = bootstrp(1000,@mean,peak(subject.Mig==1,y));
    bootPk = sort(bootPk,1);
    bootAmp = bootstrp(1000,@mean,amp(subject.Mig==1,y));
    bootAmp = sort(bootAmp,1);
    errorbar(bootPk(500),bootAmp(500),abs(diff(bootAmp([25 500]))),abs(diff(bootAmp([500 975]))),abs(diff(bootPk([25 500]))),abs(diff(bootPk([500 975]))),'o','LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','w','Color','r')
    
    figure(205)
    subplot(1,2,2)
    hold on
    bootBW = bootstrp(1000,@mean,bandwidth(subject.Mig==1,y));
    bootBW = sort(bootBW,1);
    errorbar(y,bootBW(500),abs(diff(bootBW([25 500]))),abs(diff(bootBW([500 975]))),'o','LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','w','Color','r')

end

fig = figure(205);
fig_name = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet_case/meanParamVEP_case_less28';
print(fig,fig_name,'-dpdf','-painters')

plot_meanVEP(xdata,vep(subject.Mig==0,:),'errorbars','Boot','fig_num',206,'color_mean',[0 0 0],'color_err',[0.8 0.8 0.8]);
plot_meanVEP(xdata,vep(subject.Mig==1,:),'errorbars','Boot','fig_num',206,'color_mean',[1 0 0],'color_err',[1 0.8 0.8]);

%% Plot individual VEP to check fits

figure
for i = 1:size(vep,1)
    hold on
    plot(xdata,vep(i,:),'-','Color',[0.5 0.5 0.5])
    plot(mdl_x,squeeze(sum(Gamma(i,:,:),2)),'c')
    for X = 1:nGamma
         plot(mdl_x,squeeze(Gamma(i,X,:)),['-' gammaC{X}])
         plot(squeeze(peak(i,X,:)),squeeze(amp(i,X,:)),['+' gammaC{X}])
         text(peak(i,1),amp(i,1),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[bandwidth(i,1) peak(i,1) amp(i,1)]));
         text(peak(i,2),amp(i,2),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[bandwidth(i,2) peak(i,2) amp(i,2)]));
         text(peak(i,3),amp(i,3),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[bandwidth(i,3) peak(i,3) amp(i,3)]));
         text(peak(i,4),amp(i,4),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[bandwidth(i,4) peak(i,4) amp(i,4)]));
    end
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-2 2]; ax.XLim = [0 time_end];
    xlabel(sprintf('subject %2.0f, Fit r = %2.2f',[subject.uniqueID(i) r_val(i,:)]))
    pause(1)
    clf
end


[pP100amp,tblP100amp,statsP100amp] = kruskalwallis(amp(:,2),subject.Mig);

%% plot parameters

params_tbl_predict = subjectFU;
params_tbl_predict.bw1 = bandwidth(:,1); params_tbl_predict.bw2 = bandwidth(:,2); params_tbl_predict.bw3 = bandwidth(:,3); params_tbl_predict.bw4 = bandwidth(:,4);
params_tbl_predict.pt1 = peak(:,1); params_tbl_predict.pt2 = peak(:,2); params_tbl_predict.pt3 = peak(:,3); params_tbl_predict.pt4 = peak(:,4);
params_tbl_predict.amp1 = amp(:,1); params_tbl_predict.amp2 = amp(:,2); params_tbl_predict.amp3 = amp(:,3); params_tbl_predict.amp4 = amp(:,4);

params_tbl1 = subject;
params_tbl1.bw1 = bandwidth(:,1); params_tbl1.bw2 = bandwidth(:,2); params_tbl1.bw3 = bandwidth(:,3); params_tbl1.bw4 = bandwidth(:,4);
params_tbl1.pt1 = peak(:,1); params_tbl1.pt2 = peak(:,2); params_tbl1.pt3 = peak(:,3); params_tbl1.pt4 = peak(:,4);
params_tbl1.amp1 = amp(:,1); params_tbl1.amp2 = amp(:,2); params_tbl1.amp3 = amp(:,3); params_tbl1.amp4 = amp(:,4);
