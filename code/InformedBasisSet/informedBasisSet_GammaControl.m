% Use single gamma model (can adjust to different number of gammas)


data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])
save_indivVEP_path = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/indivVEP_wFits/';
load([data_path '/neuroActive_meds.mat']) 

%% Select control subjects, first session, who do not have a history of migraine
control_subject = find(cleaned_vep_files.subjecttype=='Control' & cleaned_vep_files.med_hx___1==0);

cleaned_vep = cleaned_vep(control_subject,:);
control_vep_subjects = cleaned_vep_files(control_subject,:);
xdata = cleaned_vep{1,3}.*1000; %convert to ms

% Select first session for each participant
unique_ID = unique(control_vep_subjects.uniqueID);
control_vep = cell(length(unique_ID),3);
vep = zeros(1,length(xdata));
subject = control_vep_subjects(1,:);

counter = 1;
for x = 1:length(unique_ID)
    temp_loc = find(cell2mat(cleaned_vep(:,1))==unique_ID(x,:));
    if size(cleaned_vep{temp_loc(1),4},1)>150
        vep(counter,:) = mean(cleaned_vep{temp_loc(1),4},1).*100; %correct amplification issue with diopsys
        subject(counter,:) = control_vep_subjects(temp_loc(1),:);
        counter = counter+1;
    end
    clear temp_loc
end

%remove participants on neuroactive meds
temp = setdiff(subject.uniqueID,neuro_active_meds);
subject = subject(ismember(subject.uniqueID,temp),:);
vep = vep(ismember(subject.uniqueID,temp),:);


% shift VEP to 0 = mean response
for x = 1:size(vep)
        ydata = vep(x,:);
        vep(x,:) = ydata - mean(ydata);
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

r_val = zeros(size(control_vep,1),1);
r_val200 = zeros(size(control_vep,1),1);


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
        
figure
for i = 1:size(vep,1)
        hold on
        plot([0 500],[0 0],'--k')
        ydata = vep(i,:);
        diffY = diff([min(ydata) max(ydata)]);
        plot(xdata,ydata,'Color',[0.5 0.5 0.5])
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
        
        plot(peak75,amp75,'+b')
        
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
        plot(peak100,amp100,'+r')
       
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
        plot(peak135,amp135,'+m')
        
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
        
        plot(peak220,amp220,'+g')
        
        bw75 = 10^((80-abs(diff([peak75 peak100])))/30);

        if bw75 < 20
            bw75 = 20;
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
        
        title(sprintf('BW guess 75 = %2.2f, 100 = %2.2f, 135 = %2.2f, 220 = %2.2f',[bw75 bw100 bw135 bw220]));ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-45 45]; ax.XLim = [0 time_end];
        
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
        
        pause(1)
        clf

end

for i = 1:size(vep,1)

    ydata = vep(i,:); % averaged across trials, already corrected for diopsys amplification error above
    
    p0 = [Bw75(i,:) Peak75(i,:) Amp75(i,:) Bw100(i,:) Peak100(i,:) Amp100(i,:) Bw135(i,:) Peak135(i,:) Amp135(i,:) Bw220(i,:) Peak220(i,:) Amp220(i,:)];
    lb = [30 Peak75(i,:)-2 Amp75(i,:)*1.05 20 Peak100(i,:)-3 0.5 15 Peak135(i,:)-5 Amp135(i,:)*1.05 15 Peak220(i,:)-5 0.5]; 
    ub = [110 Peak75(i,:)+2 -0.5 110 Peak100(i,:)+3 Amp100(i,:)*1.05 110 Peak135(i,:)+5 -0.5 100 Peak220(i,:)+5 Amp220(i,:)*1.05];

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
    r = corrcoef(ydata(1,1:204),yFit(i,1:204));
    r_val200(i,:) = r(1,2);

    switch i
        case {1,9,17,25,33,41,49,57,65,73,81,89,97,105,113,121,129,137,145,153}
            fig = figure;
            fig_name = [save_indivVEP_path 'indVEP' num2str(i) 'to' num2str(i+7)];
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
    title(sprintf('%2.0f, age = %2.0f',[i subject.age_vep(i)]))
%     ylabel(sprintf('amp %2.0f, n75 %2.0f ms, p100 %2.0f ms',[peak_diff(i) lat_n75(i) lat_p100(i)]))

    subplot(4,4,j+4)
    hold on
    for X = 1:nGamma
         plot(mdl_x,gamma(X,:),['-' gammaC{X}])
    end
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-40 40]; ax.XLim = [0 time_end];
    
    switch i
        case {8,16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,136,144,152,155}
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
plot(2*ones(size(r_val)),r_val,'.','Color',[0.5 0.5 0.5])
errorbar(2,mean(r_val200),std(r_val200),'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 3];

% Plot parameters as a function of age

figure(300)
z = 1;
for x = 1:3
    for y = 1:nGamma
        subplot(3,nGamma,z)
        hold on
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    
        switch x
            case 1
                plot(subject.age_vep,bandwidth(:,y),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
                lsline
                [rr,pp] = corrcoef(subject.age_vep,bandwidth(:,y));
                rr2 = rr(1,2)^2;
                title(sprintf('bandwidth %1d',y))
                ylim([0 max(max(bandwidth))])
            case 2
                plot(subject.age_vep,peak(:,y),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
                lsline
                plot([10 20],[lb((y*3)-1) lb((y*3)-1)],'--')
                plot([10 20],[ub((y*3)-1) ub((y*3)-1)],'--')
                [rr,pp] = corrcoef(subject.age_vep,peak(:,y));
                rr2 = rr(1,2)^2;
                title(sprintf('peak time %1d',y))
                ylim([0 400])
            case 3
                plot(subject.age_vep,abs(amp(:,y)),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
                lsline
                [rr,pp] = corrcoef(subject.age_vep,amp(:,y));
                rr2 = rr(1,2)^2;
                title(sprintf('amplitude %1d',y))
                ylim([0 max(max(amp))])
        end
        xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
        z = z+1;
    end
end



%% Plot mean VEP and parameters


figure(205)
subplot(1,2,1)
hold on
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';

figure(205)
subplot(1,2,2)
hold on
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 5];

% plotting parameters

for y = 1:nGamma
     
    plot_meanVEP(mdl_x,squeeze(Gamma(:,y,:)),'errorbars','Boot','fig_num',205,'sub_plot',true,'sub_plot_num',[1 2 1],'color_mean',[0.2 0.2 0.2],'color_err',[0.8 0.8 0.8]);
    
    figure(205)
    subplot(1,2,1)
    bootPk = bootstrp(1000,@mean,peak(:,y));
    bootPk = sort(bootPk,1);
    bootAmp = bootstrp(1000,@mean,amp(:,y));
    bootAmp = sort(bootAmp,1);
    errorbar(bootPk(500),bootAmp(500),abs(diff(bootAmp([25 500]))),abs(diff(bootAmp([500 975]))),abs(diff(bootPk([25 500]))),abs(diff(bootPk([500 975]))),'o','LineWidth',2,'MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
    
    figure(205)
    subplot(1,2,2)
    hold on
    bootBW = bootstrp(1000,@mean,bandwidth(:,y));
    bootBW = sort(bootBW,1);
    errorbar(y,bootBW(500),abs(diff(bootBW([25 500]))),abs(diff(bootBW([500 975]))),'o','LineWidth',2,'MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})

end

fig = figure(205);
fig_name = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/meanParamVEP';
print(fig,fig_name,'-dpdf','-painters')

plot_meanVEP(mdl_x,squeeze(sum(Gamma,2)),'errorbars','Boot','fig_num',206,'color_mean',[0.2 0.2 0.2],'color_err',[0.8 0.8 0.8]);
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';


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
    pause
    clf
end

%% model age and VEP

mdl_age = fitlm([bandwidth peak amp],subject.age_vep);