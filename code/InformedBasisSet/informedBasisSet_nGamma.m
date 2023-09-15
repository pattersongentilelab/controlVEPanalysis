% compares gamma models with different numbers of gamma functions for vep
% analysis


data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])
saveFig_path = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/';

%% split recording into two halves and compare variance explained with a variable number of gamma functions
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
nGamma = 6;

xdata = xdata(:,1:xdata_end);
vep = vep(:,:,1:xdata_end);

meanVEP = squeeze(mean(vep,1));

r_val = zeros(size(vep,1),2);
r2 = zeros(size(vep,1),2);

%% Check for local minima and maxima for model guesses

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
Amp330 = zeros(size(vep,1),1);
Peak330 = zeros(size(vep,1),1);
Bw330 = zeros(size(vep,1),1);
Amp450 = zeros(size(vep,1),1);
Peak450 = zeros(size(vep,1),1);
Bw450 = zeros(size(vep,1),1);

for i = 1:size(vep,1)
    ydata = squeeze(vep(i,1,:));
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

    x = sum(min_loc(xdata>peak220+30 & xdata<450));
    switch x
        case 0
            amp330 = min(ydata(xdata>peak220+30 & xdata<450));
            peak330 = xdata(ydata==amp330);
            peak330 = peak330(1);
        case 1
            peak330 = min_peak(min_peak>peak220+30 & min_peak<450);
            amp330 = ydata(xdata==peak330);
        otherwise
            peak330 = min_peak(min_peak>peak220+30 & min_peak<450);
            peak330 = peak330(1);
            amp330 = ydata(xdata==peak330);
    end

    if amp330>-1
        amp330 = -1;
    end

          x = sum(max_loc(xdata>peak330+30 & xdata<500));
    switch x
        case 0
            amp450 = max(ydata(xdata>peak330+30 & xdata<500));
            peak450 = xdata(ydata==amp450);
            peak450 = peak450(1);
        case 1
            peak450 = max_peak(max_peak>peak330+30 & max_peak<500);
            amp450 = ydata(xdata==peak450);
        otherwise
            peak450 = max_peak(max_peak>peak330+30 & max_peak<500);
            peak450 = peak450(1);
            amp450 = ydata(xdata==peak450);
    end

    if amp450<1
        amp450 = 1;
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

    bw330 = 10^((250-abs(diff([peak220 peak330])))/100);
    if bw330 < 10
        bw330 = 10;
    end
    if bw330 > 80
        bw330 = 80;
    end

    bw450 = 10^((300-abs(diff([peak330 peak450])))/120);
    if bw450 < 10
        bw450 = 10;
    end
    if bw450 > 80
        bw450 = 80;
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
    Amp330(i,:) = amp330;
    Peak330(i,:) = peak330;
    Bw330(i,:) = bw330;
    Amp450(i,:) = amp450;
    Peak450(i,:) = peak450;
    Bw450(i,:) = bw450;
    
    
    figure(1)
    hold on
    plot(xdata,ydata,'Color',[0.5 0.5 0.5])
        
    plot(peak75,amp75,'+b')
    text(peak75+0.1,amp75+0.1,sprintf('bw guess = %2.2f',bw75));
    plot(peak100,amp100,'+r')
    text(peak100+0.1,amp100+0.1,sprintf('bw guess = %2.2f',bw100));
        plot(peak135,amp135,'+m')
    text(peak135+0.1,amp135+0.1,sprintf('bw guess = %2.2f',bw135));
        plot(peak220,amp220,'+g')
    text(peak220+0.1,amp220+0.1,sprintf('bw guess = %2.2f',bw220));
        plot(peak330,amp330,'+c')
    text(peak330+0.1,amp330+0.1,sprintf('bw guess = %2.2f',bw330));
        plot(peak450,amp450,'+y')
    text(peak450+0.1,amp450+0.1,sprintf('bw guess = %2.2f',bw450));

    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-45 45]; ax.XLim = [0 time_end];
    pause(1)
    clf
    
end

%% Determine fits on individual VEP data
mdl = NaN*ones(size(vep,1),nGamma,3*nGamma);
Gamma = NaN*ones(size(vep,1),nGamma,nGamma,length(mdl_x));
yFit = NaN*ones(size(vep));
bandwidth = NaN*ones(size(vep,1),nGamma,nGamma);

for I = 1:nGamma % loop across different gamma function numbers

    for i = 1:size(vep,1)

        ydata_train = squeeze(vep(i,1,:))'; % averaged across trials, already corrected for diopsys amplification error above
        ydata_test = squeeze(vep(i,2,:))'; % averaged across trials, already corrected for diopsys amplification error above
        
        p0 = [Bw75(i,:) Peak75(i,:) Amp75(i,:) Bw100(i,:) Peak100(i,:) Amp100(i,:) Bw135(i,:) Peak135(i,:) Amp135(i,:) Bw220(i,:) Peak220(i,:) Amp220(i,:) Bw330(i,:) Peak330(i,:) Amp330(i,:) Bw450(i,:) Peak450(i,:) Amp450(i,:)];
        lb = [max([30 Bw75(i,:)-5]) Peak75(i,:)-2 Amp75(i,:)*1.1 max([20 Bw100(i,:)-5]) Peak100(i,:)-3 Amp100(i,:)*0.9 max([15 Bw135(i,:)-5]) Peak135(i,:)-5 Amp135(i,:)*1.1 max([15 Bw220(i,:)-5]) Peak220(i,:)-5 Amp220(i,:)*0.9 max([15 Bw330(i,:)-5]) Peak330(i,:)-5 Amp330(i,:)*1.1 max([15 Bw450(i,:)-5]) Peak450(i,:)-5 Amp450(i,:)*0.9]; 
        ub = [min([110 Bw75(i,:)+5]) Peak75(i,:)+2 Amp75(i,:)*0.9 min([110 Bw100(i,:)+5]) Peak100(i,:)+3 Amp100(i,:)*1.1 min([110 Bw135(i,:)+5]) Peak135(i,:)+5 -Amp135(i,:)*0.9 min([100 Bw220(i,:)+5]) Peak220(i,:)+5 Amp220(i,:)*1.1 min([100 Bw330(i,:)+5]) Peak330(i,:)+5 -Amp330(i,:)*0.9 min([100 Bw450(i,:)+5]) Peak450(i,:)+5 Amp450(i,:)*1.1];
        
        p0 = p0(1:I*3);
        lb = lb(1:I*3);
        ub = ub(1:I*3);

        myFx = @(p) sqrt(sum((ydata_train - gammaVEP_model(xdata,p,I)).^2));
        mdl(i,I,1:length(p0)) = fmincon(myFx,p0,[],[],[],[],lb,ub);
        [vep_fit,gamma] = gammaVEP_model(mdl_x,squeeze(mdl(i,I,:)),I);

        for z = 1:I
            bandwidth(i,I,z) = gamma_bandwidth(mdl_x,gamma(z,:));
            Gamma(i,I,z,:) = gamma(z,:);
        end

        [yFit(i,I,:)] = gammaVEP_model(xdata,squeeze(mdl(i,I,:)),I);
        r = corrcoef(ydata_test,squeeze(yFit(i,I,:)));
        r_val(i,I) = r(1,2);
        r2(i,I) = r(1,2)^2;
    end
    % Plot R values
    figure(202)
    hold on
    r_val95 = bootstrp(1000,@median,r_val(:,I));
    plot(I*ones(size(r_val(:,I))),r_val(:,I),'.','Color',gammaC{I})
    errorbar(I,prctile(r_val95,50),abs(diff(prctile(r_val95,[2.5 50]))),abs(diff(prctile(r_val95,[50 97.5]))),'o','LineWidth',2,'MarkerFaceColor',gammaC{I},'MarkerEdgeColor','w','Color',gammaC{I})
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 nGamma+1];
end

fig = figure(202);
fig_name = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/meanFitQual_nGamma';
print(fig,fig_name,'-dpdf','-painters')

peak = mdl(:,:,2:3:end);
amp = mdl(:,:,3:3:end);



% plotting parameters for the full 6 component gamma model

for y = 1:nGamma
     
    plot_meanVEP(mdl_x,squeeze(Gamma(:,nGamma,y,:)),'errorbars','Boot','fig_num',205,'sub_plot',true,'sub_plot_num',[1 3 1],'color_mean',[0.2 0.2 0.2],'color_err',[0.8 0.8 0.8]);
    
    pk = squeeze(peak(:,nGamma,:));
    am = squeeze(amp(:,nGamma,:));
    bw = squeeze(bandwidth(:,nGamma,:));
    
    figure(205)
    subplot(1,3,1)
    bootPk = bootstrp(1000,@mean,pk(:,y));
    bootPk = sort(bootPk,1);
    bootAmp = bootstrp(1000,@mean,am(:,y));
    bootAmp = sort(bootAmp,1);
    errorbar(bootPk(500),bootAmp(500),abs(diff(bootAmp([25 500]))),abs(diff(bootAmp([500 975]))),abs(diff(bootPk([25 500]))),abs(diff(bootPk([500 975]))),'o','LineWidth',2,'MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})
    
    figure(205)
    subplot(1,3,2)
    hold on
    bootBW = bootstrp(1000,@mean,bw(:,y));
    bootBW = sort(bootBW,1);
    errorbar(y,bootBW(500),abs(diff(bootBW([25 500]))),abs(diff(bootBW([500 975]))),'o','LineWidth',2,'MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})

end

plot_meanVEP(mdl_x,squeeze(sum(squeeze(Gamma(:,end,:,:)),2)),'errorbars','Boot','fig_num',205,'sub_plot',true,'sub_plot_num',[1 3 3],'color_mean',[0.2 0.2 0.2],'color_err',[0.8 0.8 0.8]);
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';

fig = figure(205);
fig_name = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/meanParamVEPnGamma';
print(fig,fig_name,'-dpdf','-painters')