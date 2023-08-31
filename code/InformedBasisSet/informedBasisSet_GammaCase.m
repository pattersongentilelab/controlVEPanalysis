% Use single gamma model (can adjust to different number of gammas)

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])

%% select case subjects who have VEP <28 days and >28 days post injury

case_subject = find(cleaned_vep_files.subjecttype=='Case');
cleaned_vep = cleaned_vep(case_subject,:);
case_vep_subjects = cleaned_vep_files(case_subject,:);
xdata = cleaned_vep{1,3}.*1000; %convert to ms

unique_ID = unique(case_vep_subjects.uniqueID);
case_vep = cell(length(unique_ID),3);
vep = zeros(length(unique_ID),length(xdata));
subject = control_vep_subjects(1,:);

for x = 1:length(unique_ID)
    temp_loc = find(cell2mat(cleaned_vep(:,1))==unique_ID(x,:));
    case_vep(x,:) = cleaned_vep(temp_loc(1),[1 2 4]);
    vep(x,:) = mean(cleaned_vep{temp_loc(1),4},1).*100; %correct amplification issue with diopsys
    subject(x,:) = control_vep_subjects(temp_loc(1),:);
    clear temp_loc
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
bootstat = bootstrp(1000,@mean,vep); boostat = sort(bootstat,1);
bootVEP = bootstat([25 500 975],:);

r_val = zeros(size(control_vep,1),1);
r2 = zeros(size(control_vep,1),1);


 % Determine fit on mean VEP data across individuals
ydata = meanVEP;
% Model guess for average model
p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 25 300 max(ydata)];
lb = [10 65 min(ydata)*1.2 10 70 0.5 10 110 min(ydata)*1.2 10 200 0.5]; 
ub = [60 90 -0.5 60 120 max(ydata)*1.2 60 160 -0.5 60 400 max(ydata)*1.2];

myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,nGamma)).^2));
Mdl = fmincon(myFx,p0,[],[],[],[],lb,ub);
[yFit_m,gammaM] = gammaVEP_model(xdata,Mdl,nGamma);

figure(205)
hold on
errorbar(xdata,bootVEP(2,:),abs(diff(bootVEP(1:2,:),1)),abs(diff(bootVEP(2:3,:),1)),'.k')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
plot(xdata,yFit_m,['-' gammaC{nGamma}],'LineWidth',2)


%% Determine fits on individual VEP data
mdl = zeros(size(control_vep,1),3*nGamma);
gamma = zeros(size(control_vep,1),length(xdata));
yFit = zeros(size(vep));
bandwidth = NaN*ones(size(vep,1),nGamma);
    
figure
for i = 1:size(control_vep,1)

    ydata = vep(i,:); % averaged across trials, already corrected for diopsys amplification error above

    p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 25 300 max(ydata)];
    lb = [10 65 min(ydata)*1.2 10 70 0.5 10 110 min(ydata)*1.2 10 200 0.5]; 
    ub = [60 90 -0.5 60 120 max(ydata)*1.2 60 160 -0.5 60 400 max(ydata)*1.2];

    myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,nGamma)).^2));
    mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [vep_fit,gamma] = gammaVEP_model(mdl_x,mdl(i,:),nGamma);

    zz = 1;
    for z = 1:nGamma
        bandwidth(i,z,:) = gamma_bandwidth(0,500,mdl(i,zz:zz+1));
        zz = zz+3;
    end

    [yFit(i,:)] = gammaVEP_model(xdata,mdl(i,:),nGamma);
    r = corrcoef(ydata,yFit(i,:));
    r_val(i,:) = r(1,2);
    r2(i,:) = r(1,2)^2;

    switch i
        case {1,9,17,25,33,41,49,57,65,73,81,89,97,105,113,121,129,137,145,153}
            figure
            j = 1;
        case {5,13,21,29,37,45,53,61,69,77,85,93,101,109,117,125,133,141,149}
            j = 9;
        otherwise
            j = j+1;
    end

    subplot(4,4,j)
    plot(xdata,vep(i,:),'.k')
    hold on
    plot(xdata,yFit(i,:),['-' gammaC{nGamma}],'LineWidth',2)
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    xlabel(sprintf('r = %2.2f',r_val(i)))
    title(sprintf('%2.0f, age = %2.0f',[i subject.age_vep(i)]))

    subplot(4,4,j+4)
    hold on
    for X = 1:nGamma
         plot(mdl_x,gamma(X,:),['-' gammaC{X}])
    end
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
end

peak = mdl(:,2:3:end);
amp = mdl(:,3:3:end);   

% plotting parameters
figure(200)
for x = 1:3
    subplot(1,3,x)
    hold on
    for y = 1:nGamma
    
        switch x
            case 1
                bootstat = bootstrp(1000,@mean,bandwidth(:,y));
                title('bandwidth')
            case 2
                bootstat = bootstrp(1000,@mean,peak(:,y));                
                title('peak time')
            case 3
                bootstat = bootstrp(1000,@mean,amp(:,y));
                title('amplitude')
        end
        boostat = sort(bootstat,1);
        boot = bootstat([25 500 975]);
        errorbar(y,boot(2),abs(diff(boot(1:2))),abs(diff(boot(2:3))),'o','LineWidth',2,'MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
        
    end
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:nGamma; ax.XLim = [0 nGamma+1];
end

    
% Plot R values
figure(202)
hold on
plot(ones(size(r_val)),r_val,'.','Color',[0.5 0.5 0.5])
errorbar(1,mean(r_val),std(r_val),'o','LineWidth',2,'MarkerFaceColor',gammaC{nGamma},'MarkerEdgeColor',gammaC{nGamma},'Color',gammaC{nGamma})
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 2];

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
                ylim([min(min(amp)) max(max(amp))])
        end
        xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
        z = z+1;
    end
end

% models for age and VEP
% X = subject.age_vep;
% Y = [bandwidth peak amp];
% 
% mdl_gamma = fitlm(X,Y);

%% Compare to diopsys peak analysis
lat_n75 = subject.Left_Cursor_Lat(subject.Left_Cursor_Lat~=0);
lat_p100 = subject.Right_Cursor_Lat(subject.Left_Cursor_Lat~=0);

peak_n75 = subject.Left_Cursor_Amp(subject.Left_Cursor_Lat~=0);
peak_p100 = subject.Right_Cursor_Amp(subject.Left_Cursor_Lat~=0);

peak_diff = abs(diff([peak_n75 peak_p100],[],2));

peak_diffIBS = abs(diff(amp(subject.Left_Cursor_Lat~=0,1:2),[],2));

figure(400)
% Peak-to-peak amplitude
subplot(2,2,[1 2])
hold on
plot(peak_diff,peak_diffIBS,'.k')
lsline
[rr,pp] = corrcoef(peak_diff,peak_diffIBS);
title(sprintf('peak-to-peak amplitude, r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
xlabel('ISCEV peak analysis')
ylabel('informed basis set')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off';

% peak time
subplot(2,2,3)
hold on
plot(lat_n75,peak(subject.Left_Cursor_Lat~=0,1),'.k')
lsline
[rr,pp] = corrcoef(lat_n75,peak(subject.Left_Cursor_Lat~=0,1));
title(sprintf('N75 peak time, r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
xlabel('ISCEV peak analysis')
ylabel('informed basis set')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
ax.XLim = [0 150];  ax.YLim = [0 150];

subplot (2,2,4)
hold on
plot(lat_p100,peak(subject.Left_Cursor_Lat~=0,2),'.k')
lsline
[rr,pp] = corrcoef(lat_p100,peak(subject.Left_Cursor_Lat~=0,2));
title(sprintf('P100 peak time, r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
xlabel('ISCEV peak analysis')
ylabel('informed basis set')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
ax.XLim = [0 200];  ax.YLim = [0 200];