% Use single gamma model (can adjust to different number of gammas)


data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])
save_indivVEP_path = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/indivVEP_wFits/';
load([data_path '/neuroActive_meds.mat']) 

%% Select control subjects, first session
control_subject = find(cleaned_vep_files.subjecttype=='Control');

cleaned_vep = cleaned_vep(control_subject,:);
control_vep_subjects = cleaned_vep_files(control_subject,:);
xdata = cleaned_vep{1,3}.*1000; %convert to ms

% Select first session for each participant
vep = zeros(1,length(xdata));
subject = control_vep_subjects(1,:);
uniqueID = unique(control_vep_subjects.uniqueID);

counter = 1;
for x = 1:length(uniqueID)
    temp_loc = find(cell2mat(cleaned_vep(:,1))==uniqueID(x));
    if size(cleaned_vep{temp_loc(1),4},1)>150
        vep(counter,:) = mean(cleaned_vep{temp_loc(1),4},1).*100; %correct amplification issue with diopsys
        subject(counter,:) = control_vep_subjects(temp_loc(1),:);
        counter = counter+1;
    end
    clear temp_loc
end
clear uniqueID



% shift VEP to 0 = mean response
for x = 1:size(vep)
        ydata = vep(x,:);
        ydata = ydata - mean(ydata);
%         ydata = ydata./max(abs(ydata)); % normalize by absolute maximum response
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
        if amp75>-1
            amp75 = -1;
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
        
        if amp100<1
            amp100 = 1;
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
        
        if amp135>-1
            amp135 = -1;
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
 
        if amp220<1
            amp220 = 1;
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

    ydata = vep(i,:); % averaged across trials, already corrected for diopsys amplification error above
    
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
        case {1,9,17,25,33,41,49,57,65,73,81,89,97,105,113,121,129,137,145}
            fig = figure;
            j = 1;
        case {5,13,21,29,37,45,53,61,69,77,85,93,101,109,117,125,133,141,149}
            j = 9;
        otherwise
            j = j+1;
    end

    subplot(4,4,j)
    plot(xdata,vep(i,:),'-','Color',[0.5 0.5 0.5])
    hold on
    plot(xdata,yFit(i,:),'-k','LineWidth',2)
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 time_end]; ax.YLim = [-40 40];
    xlabel(sprintf('r = %2.2f',r_val(i)))
    title(sprintf('%2.0f, C%2.0f, age = %2.0f',[i subject.uniqueID(i) subject.age_vep(i)]))
    subplot(4,4,j+4)
    hold on
    for X = 1:nGamma
         plot(mdl_x,gamma(X,:),['-' gammaC{X}])
    end
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-40 40]; ax.XLim = [0 time_end];
    
    switch i
        case {8,16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,136,144,151}
        fig_name = [save_indivVEP_path 'indVEP' num2str(i-7) 'to' num2str(i)];
        print(fig,fig_name,'-dpdf','-painters')
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
subplot(1,3,1)
hold on
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';

figure(205)
subplot(1,3,2)
hold on
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 5];

% plotting parameters

for y = 1:nGamma
     
    plot_meanVEP(mdl_x,squeeze(Gamma(:,y,:)),'errorbars','Boot','fig_num',205,'sub_plot',true,'sub_plot_num',[1 3 1],'color_mean',[0.2 0.2 0.2],'color_err',[0.8 0.8 0.8]);
    
    figure(205)
    subplot(1,3,1)
    bootPk = bootstrp(1000,@mean,peak(:,y));
    bootPk = sort(bootPk,1);
    bootAmp = bootstrp(1000,@mean,amp(:,y));
    bootAmp = sort(bootAmp,1);
    errorbar(bootPk(500),bootAmp(500),abs(diff(bootAmp([25 500]))),abs(diff(bootAmp([500 975]))),abs(diff(bootPk([25 500]))),abs(diff(bootPk([500 975]))),'o','LineWidth',2,'MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})
    
    figure(205)
    subplot(1,3,2)
    hold on
    bootBW = bootstrp(1000,@mean,bandwidth(:,y));
    bootBW = sort(bootBW,1);
    errorbar(y,bootBW(500),abs(diff(bootBW([25 500]))),abs(diff(bootBW([500 975]))),'o','LineWidth',2,'MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})

end

plot_meanVEP(mdl_x,squeeze(sum(Gamma,2)),'errorbars','Boot','fig_num',205,'sub_plot',true,'sub_plot_num',[1 3 3],'color_mean',[0.2 0.2 0.2],'color_err',[0.8 0.8 0.8]);
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';

fig = figure(205);
fig_name = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/meanParamVEP';
print(fig,fig_name,'-dpdf','-painters')




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
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-40 40]; ax.XLim = [0 time_end];
    xlabel(sprintf('subject %2.0f, Fit r = %2.2f',[subject.uniqueID(i) r_val(i,:)]))
    pause(1)
    clf
end

%% comparison of peak analysis (diopsys) to difference-of-gammas model

% Compile diopsys measurements
load([data_path '/VEPsubject_masterList111220.mat'])

diopsys = subject(:,1);
diopsys.n75_pt = zeros(height(diopsys),1);
diopsys.p100_pt = zeros(height(diopsys),1);
diopsys.delta_amp = zeros(height(diopsys),1);
diopsys.sessNum = zeros(height(diopsys),1);

for x = 1:height(diopsys)
    temp = VEPsubject_data(VEPsubject_data.uniqueID==diopsys.uniqueID(x) & VEPsubject_data.REDCap_testdate==subject.REDCap_testdate(x),:);
    temp = temp(temp.Delta_Amp~=0 & temp.Left_Cursor_Lat~=0 & temp.Right_Cursor_Lat~=0,:);
    diopsys.n75_pt(x) = mean(temp.Left_Cursor_Lat);
    diopsys.p100_pt(x) = mean(temp.Right_Cursor_Lat);
    diopsys.delta_amp(x) = mean(temp.Delta_Amp);
    diopsys.sessNum(x) = height(temp);
end


n75_pt = diopsys.n75_pt;
p100_pt = diopsys.p100_pt;
n75p100_amp = diopsys.delta_amp;

n75_ptG = peak(:,1);
p100_ptG = peak(:,2);
n75p100_ampG = abs(diff(amp(:,1:2),1,2));

% Bland Altman Plots

fig = figure(302);
subplot(1,3,1)
diff_n75pt = diff([n75_pt n75_ptG],[],2);
mean_n75pt = mean([n75_pt n75_ptG],2);
hold on
plot(mean_n75pt,diff_n75pt,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','Color',[0.5 0.5 0.5])
plot([min(mean_n75pt) max(mean_n75pt)],[mean(diff_n75pt) mean(diff_n75pt)],'--k')
plot([min(mean_n75pt) max(mean_n75pt)],[mean(diff_n75pt)+std(diff_n75pt)*1.96 mean(diff_n75pt)+std(diff_n75pt)*1.96],'--k')
plot([min(mean_n75pt) max(mean_n75pt)],[mean(diff_n75pt)+std(diff_n75pt)*-1.96 mean(diff_n75pt)+std(diff_n75pt)*-1.96],'--k')
text(max(mean_n75pt)-2,mean(diff_n75pt)+std(diff_n75pt)*-1.96-2,sprintf('%2.1f',mean(diff_n75pt)+std(diff_n75pt)*-1.96));
text(max(mean_n75pt)-2,mean(diff_n75pt)+std(diff_n75pt)*1.96+2,sprintf('%2.1f',mean(diff_n75pt)+std(diff_n75pt)*1.96));
title('N75 peak time')
xlabel('Average of methods')
ylabel('Difference between methods')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75pt), max(mean_n75pt)]; ax.YLim = [-50,50];


subplot(1,3,2)
diff_p100pt = diff([p100_pt p100_ptG],[],2);
mean_p100pt = mean([p100_pt p100_ptG],2);
hold on
plot(mean_p100pt,diff_p100pt,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','Color',[0.5 0.5 0.5])
plot([min(mean_p100pt) max(mean_p100pt)],[mean(diff_p100pt) mean(diff_p100pt)],'--k')
plot([min(mean_p100pt) max(mean_p100pt)],[mean(diff_p100pt)+std(diff_p100pt)*1.96 mean(diff_p100pt)+std(diff_p100pt)*1.96],'--k')
plot([min(mean_p100pt) max(mean_p100pt)],[mean(diff_p100pt)+std(diff_p100pt)*-1.96 mean(diff_p100pt)+std(diff_p100pt)*-1.96],'--k')
text(max(mean_p100pt)-2,mean(diff_p100pt)+std(diff_p100pt)*-1.96-2,sprintf('%2.1f',mean(diff_p100pt)+std(diff_p100pt)*-1.96));
text(max(mean_p100pt)-2,mean(diff_p100pt)+std(diff_p100pt)*1.96+2,sprintf('%2.1f',mean(diff_p100pt)+std(diff_p100pt)*1.96));
title('P100 peak time')
xlabel('Average of methods')
ylabel('Difference between methods')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_p100pt), max(mean_p100pt)]; ax.YLim = [-50,50];

subplot(1,3,3)
diff_amp = diff([n75p100_amp n75p100_ampG],[],2);
mean_amp = mean([n75p100_amp n75p100_ampG],2);
hold on
plot(mean_amp,diff_amp,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','Color',[0.5 0.5 0.5])
plot([min(mean_amp) max(mean_amp)],[mean(diff_amp) mean(diff_amp)],'--k')
plot([min(mean_amp) max(mean_amp)],[mean(diff_amp)+std(diff_amp)*1.96 mean(diff_amp)+std(diff_amp)*1.96],'--k')
plot([min(mean_amp) max(mean_amp)],[mean(diff_amp)+std(diff_amp)*-1.96 mean(diff_amp)+std(diff_amp)*-1.96],'--k')
text(max(mean_amp)-2,mean(diff_amp)+std(diff_amp)*-1.96-2,sprintf('%2.1f',mean(diff_amp)+std(diff_amp)*-1.96));
text(max(mean_amp)-2,mean(diff_amp)+std(diff_amp)*1.96+2,sprintf('%2.1f',mean(diff_amp)+std(diff_amp)*1.96));
title('N75-P100 amplitude')
xlabel('Average of methods')
ylabel('Difference between methods')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_amp), max(mean_amp)]; ax.YLim = [-50,50];

fig_name = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/DoG_vs_PeakAnalysis';
print(fig,fig_name,'-dpdf','-painters')

%% age and VEP



params_tbl = table([subject.age_vep],'VariableNames',{'age'});
params_tbl.sex = subject.sex_master-1;
params_tbl.concHx = subject.conc_hx_yn;
params_tbl.bw1 = bandwidth(:,1); params_tbl.bw2 = bandwidth(:,2); params_tbl.bw3 = bandwidth(:,3); params_tbl.bw4 = bandwidth(:,4);
params_tbl.pt1 = peak(:,1); params_tbl.pt2 = peak(:,2); params_tbl.pt3 = peak(:,3); params_tbl.pt4 = peak(:,4);
params_tbl.amp1 = amp(:,1); params_tbl.amp2 = amp(:,2); params_tbl.amp3 = amp(:,3); params_tbl.amp4 = amp(:,4);

% Plot parameters as a function of age, no significant past medical history
% only

figure(300)

z = 1;
for x = 1:(3*nGamma)
    subplot(nGamma,3,z)
    hold on
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';

    switch x
        case 1
            plot(params_tbl.age,params_tbl.bw1,'o','MarkerFaceColor',gammaC{1},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.bw1);
            rr2 = rr(1,2)^2;
            title('bandwidth 1')
            ylim([0 max(params_tbl.bw4)])
        case 2
            plot(params_tbl.age,params_tbl.pt1,'o','MarkerFaceColor',gammaC{1},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.pt1);
            rr2 = rr(1,2)^2;
            title('peak time 1')
            ylim([0 400])
        case 3
            plot(params_tbl.age,params_tbl.amp1,'o','MarkerFaceColor',gammaC{1},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.amp1);
            rr2 = rr(1,2)^2;
            title('amplitude 1')
            ylim([min(params_tbl.amp1) 0])
        case 4
            plot(params_tbl.age,params_tbl.bw2,'o','MarkerFaceColor',gammaC{2},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.bw2);
            rr2 = rr(1,2)^2;
            title('bandwidth 2')
            ylim([0 max(params_tbl.bw4)])
        case 5
            plot(params_tbl.age,params_tbl.pt2,'o','MarkerFaceColor',gammaC{2},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.pt2);
            rr2 = rr(1,2)^2;
            title('peak time 2')
            ylim([0 400])
        case 6
            plot(params_tbl.age,params_tbl.amp2,'o','MarkerFaceColor',gammaC{2},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.amp2);
            rr2 = rr(1,2)^2;
            title('amplitude 2')
            ylim([0 max(params_tbl.amp2)])
        case 7
            plot(params_tbl.age,params_tbl.bw3,'o','MarkerFaceColor',gammaC{3},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.bw3);
            rr2 = rr(1,2)^2;
            title('bandwidth 3')
            ylim([0 max(params_tbl.bw4)])
        case 8
            plot(params_tbl.age,params_tbl.pt3,'o','MarkerFaceColor',gammaC{3},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.pt3);
            rr2 = rr(1,2)^2;
            title('peak time 3')
            ylim([0 400])
        case 9
            plot(params_tbl.age,params_tbl.amp3,'o','MarkerFaceColor',gammaC{3},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.amp3);
            rr2 = rr(1,2)^2;
            ylim([min(params_tbl.amp3) 0])
            title('amplitude 3')
        case 10
            plot(params_tbl.age,params_tbl.bw4,'o','MarkerFaceColor',gammaC{4},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.bw4);
            rr2 = rr(1,2)^2;
            title('bandwidth 4')
            ylim([0 max(params_tbl.bw4)])
        case 11
            plot(params_tbl.age,params_tbl.pt4,'o','MarkerFaceColor',gammaC{4},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.pt4);
            rr2 = rr(1,2)^2;
            title('peak time 4')
            ylim([0 400])
        case 12
            plot(params_tbl.age,params_tbl.amp4,'o','MarkerFaceColor',gammaC{4},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(params_tbl.age,params_tbl.amp4);
            rr2 = rr(1,2)^2;
            title('amplitude 4')
            ylim([0 max(params_tbl.amp4)])

    end
    xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
    z = z+1;
end

% find patients with no chronic medical conditions and not on neuro-active
% medications
noMedHx = params_tbl(~ismember(subject.uniqueID,neuro_active_meds) & subject.med_hx___1~=1 & subject.med_hx___2~=1 ...
& subject.med_hx___3~=1 & subject.med_hx___4~=1 & subject.med_hx___6~=1 & subject.med_hx___15~=1 & subject.med_hx___16~=1 & subject.med_hx___17==0 ...
& subject.med_hx___18~=1 & subject.med_hx___19~=1 & subject.med_hx___20~=1 & subject.med_hx___21~=1 & subject.med_hx___23~=1 & subject.med_hx___24~=1 ...
& subject.med_hx___25~=1 & subject.med_hx___26~=1 & subject.med_hx___27~=1 & subject.med_hx___28~=1,:);

vep_noMedHx = vep(~ismember(subject.uniqueID,neuro_active_meds) & subject.med_hx___1~=1 & subject.med_hx___2~=1 ...
& subject.med_hx___3~=1 & subject.med_hx___4~=1 & subject.med_hx___6~=1 & subject.med_hx___15~=1 & subject.med_hx___16~=1 & subject.med_hx___17==0 ...
& subject.med_hx___18~=1 & subject.med_hx___19~=1 & subject.med_hx___20~=1 & subject.med_hx___21~=1 & subject.med_hx___23~=1 & subject.med_hx___24~=1 ...
& subject.med_hx___25~=1 & subject.med_hx___26~=1 & subject.med_hx___27~=1 & subject.med_hx___28~=1,:);

mdl_age_noMedHx = fitlm(noMedHx,'age ~ bw1 + bw2 + bw3 + bw4 + pt1 + pt2 + pt3 + pt4 + amp1 + amp2 + amp3 + amp4');

mdl_sex_noMedHx = fitglm(noMedHx,'sex ~ bw1 + bw2 + bw3 + bw4 + pt1 + pt2 + pt3 + pt4 + amp1 + amp2 + amp3 + amp4','Distribution','binomial');

mdl_concHx_noMedHx = fitglm(noMedHx,'concHx ~ bw1 + bw2 + bw3 + bw4 + pt1 + pt2 + pt3 + pt4 + amp1 + amp2 + amp3 + amp4','Distribution','binomial');
% also remove those with remote history of concussion
% find patients with no chronic medical conditions and not on neuro-active
% medications
noMedHxConc = params_tbl(~ismember(subject.uniqueID,neuro_active_meds) & subject.med_hx___1~=1 & subject.med_hx___2~=1 ...
& subject.med_hx___3~=1 & subject.med_hx___4~=1 & subject.med_hx___6~=1 & subject.med_hx___15~=1 & subject.med_hx___16~=1 & subject.med_hx___17==0 ...
& subject.med_hx___18~=1 & subject.med_hx___19~=1 & subject.med_hx___20~=1 & subject.med_hx___21~=1 & subject.med_hx___23~=1 & subject.med_hx___24~=1 ...
& subject.med_hx___25~=1 & subject.med_hx___26~=1 & subject.med_hx___27~=1 & subject.med_hx___28~=1 & subject.conc_hx_yn~=1,:);

vep_noMedHxConc = vep(~ismember(subject.uniqueID,neuro_active_meds) & subject.med_hx___1~=1 & subject.med_hx___2~=1 ...
& subject.med_hx___3~=1 & subject.med_hx___4~=1 & subject.med_hx___6~=1 & subject.med_hx___15~=1 & subject.med_hx___16~=1 & subject.med_hx___17==0 ...
& subject.med_hx___18~=1 & subject.med_hx___19~=1 & subject.med_hx___20~=1 & subject.med_hx___21~=1 & subject.med_hx___23~=1 & subject.med_hx___24~=1 ...
& subject.med_hx___25~=1 & subject.med_hx___26~=1 & subject.med_hx___27~=1 & subject.med_hx___28~=1 & subject.conc_hx_yn~=1,:);



% Plot parameters as a function of age, no significant past medical history
% only

figure(301)

z = 1;
for x = 1:(3*nGamma)
    subplot(nGamma,3,z)
    hold on
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';

    switch x
        case 1
            plot(noMedHxConc.age,noMedHxConc.bw1,'o','MarkerFaceColor',gammaC{1},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.bw1);
            rr2 = rr(1,2)^2;
            title('bandwidth 1')
            ylim([0 max(noMedHxConc.bw4)])
        case 2
            plot(noMedHxConc.age,noMedHxConc.pt1,'o','MarkerFaceColor',gammaC{1},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.pt1);
            rr2 = rr(1,2)^2;
            title('peak time 1')
            ylim([0 400])
        case 3
            plot(noMedHxConc.age,noMedHxConc.amp1,'o','MarkerFaceColor',gammaC{1},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.amp1);
            rr2 = rr(1,2)^2;
            title('amplitude 1')
            ylim([min(noMedHxConc.amp1) 0])
        case 4
            plot(noMedHxConc.age,noMedHxConc.bw2,'o','MarkerFaceColor',gammaC{2},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.bw2);
            rr2 = rr(1,2)^2;
            title('bandwidth 2')
            ylim([0 max(noMedHxConc.bw4)])
        case 5
            plot(noMedHxConc.age,noMedHxConc.pt2,'o','MarkerFaceColor',gammaC{2},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.pt2);
            rr2 = rr(1,2)^2;
            title('peak time 2')
            ylim([0 400])
        case 6
            plot(noMedHxConc.age,noMedHxConc.amp2,'o','MarkerFaceColor',gammaC{2},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.amp2);
            rr2 = rr(1,2)^2;
            title('amplitude 2')
            ylim([0 max(noMedHxConc.amp2)])
        case 7
            plot(noMedHxConc.age,noMedHxConc.bw3,'o','MarkerFaceColor',gammaC{3},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.bw3);
            rr2 = rr(1,2)^2;
            title('bandwidth 3')
            ylim([0 max(noMedHxConc.bw4)])
        case 8
            plot(noMedHxConc.age,noMedHxConc.pt3,'o','MarkerFaceColor',gammaC{3},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.pt3);
            rr2 = rr(1,2)^2;
            title('peak time 3')
            ylim([0 400])
        case 9
            plot(noMedHxConc.age,noMedHxConc.amp3,'o','MarkerFaceColor',gammaC{3},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.amp3);
            rr2 = rr(1,2)^2;
            ylim([min(noMedHxConc.amp3) 0])
            title('amplitude 3')
        case 10
            plot(noMedHxConc.age,noMedHxConc.bw4,'o','MarkerFaceColor',gammaC{4},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.bw4);
            rr2 = rr(1,2)^2;
            title('bandwidth 4')
            ylim([0 max(noMedHxConc.bw4)])
        case 11
            plot(noMedHxConc.age,noMedHxConc.pt4,'o','MarkerFaceColor',gammaC{4},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.pt4);
            rr2 = rr(1,2)^2;
            title('peak time 4')
            ylim([0 400])
        case 12
            plot(noMedHxConc.age,noMedHxConc.amp4,'o','MarkerFaceColor',gammaC{4},'MarkerEdgeColor','w','Color',gammaC{y})
            lsline
            [rr,pp] = corrcoef(noMedHxConc.age,noMedHxConc.amp4);
            rr2 = rr(1,2)^2;
            title('amplitude 4')
            ylim([0 max(noMedHxConc.amp4)])

    end
    xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
    z = z+1;
end


params = [params_tbl.bw1 params_tbl.bw2 params_tbl.bw3 params_tbl.bw4 params_tbl.pt1 params_tbl.pt2 params_tbl.pt3 params_tbl.pt4 abs(params_tbl.amp1) params_tbl.amp2 abs(params_tbl.amp3) params_tbl.amp4];
xParam = corrcoef(params);
