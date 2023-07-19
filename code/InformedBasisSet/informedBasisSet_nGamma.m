% compares gamma models with different numbers of gamma functions for vep
% analysis

%% Determine amount of variance explained on half of the trials, predicted by the first half

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])

% select control subjects
control_subject = find(cleaned_vep_files.subjecttype=='Control');
cleaned_vep = cleaned_vep(control_subject,:);
control_vep_subjects = cleaned_vep_files(control_subject,:);

% Select first session for each participant
xdata = cleaned_vep{1,3}.*1000; %convert to ms
unique_ID = unique(control_vep_subjects.uniqueID);
control_vep = cell(length(unique_ID),3);
vep = zeros(length(unique_ID),length(xdata));

for x = 1:length(unique_ID)
    temp_loc = find(cell2mat(cleaned_vep(:,1))==unique_ID(x,:));
    control_vep(x,:) = cleaned_vep(temp_loc(1),[1 2 4]);
    vep(x,:) = mean(cleaned_vep{temp_loc(1),4},1).*100; %correct amplification issue with diopsys
    clear temp_loc
end

meanVEP = mean(vep,1);

gammaC = {'b','r','k','g','m','c','y'};
nGamma = 7;
mdl_x = 0:0.1:500;
r2 = zeros(6,size(control_vep,1));
trialR2 = zeros(size(control_vep,1),nGamma,100);
Peak = cell(nGamma);
Amp = cell(nGamma);
Bandwidth = cell(nGamma);

figure(205)
plot(xdata,meanVEP,'.k')
hold on
ax=gca; ax.TickDir = 'out'; ax.Box = 'off';

% compare fits for different gamma models
for I = 1:nGamma
    
    mdl = zeros(size(control_vep,1),3*I);
    gamma = zeros(I,size(control_vep,1),length(xdata));
    yFit = zeros(size(vep));
    bandwidth = NaN*ones(size(vep,1),I);

    % Determine fits on individual VEP data
    
    figure
    for i = 1:size(control_vep,1)
        
        Ydata = control_vep{i,3}.*100; % with separate trials, correct amplification issue with diopsys
        ydata = vep(i,:); % averaged across trials, already corrected for diopsys amplification above
        
        switch I
            case 1
                p0 = [35 75 min(ydata)];
                lb = [0 0 -Inf]; 
                ub = [Inf Inf 0];
            case 2
                p0 = [35 75 min(ydata) 25 100 max(ydata)];
                lb = [0 0 -Inf 0 0 0]; 
                ub = [Inf Inf 0 Inf Inf Inf];
            case 3
                p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata)];
                lb = [0 0 -Inf 0 0 0 0 0 -Inf]; 
                ub = [Inf Inf 0 Inf Inf Inf Inf Inf 0];
            case 4
                p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata)];
                lb = [0 0 -Inf 0 0 0 0 0 -Inf 0 0 0]; 
                ub = [Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf];
            case 5
                p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata)];
                lb = [0 0 -Inf 0 0 0 0 0 -Inf 0 0 0 0 0 -Inf]; 
                ub = [Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf Inf Inf 0];
            case 6
                p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata) 30 400 max(ydata)];
                lb = [0 0 -Inf 0 0 0 0 0 -Inf 0 0 0 0 0 -Inf 0 0 0]; 
                ub = [Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf];
            case 7
                p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata) 30 400 max(ydata) 30 400 min(ydata)];
                lb = [0 0 -Inf 0 0 0 0 0 -Inf 0 0 0 0 0 -Inf 0 0 0 0 0 0 -Inf]; 
                ub = [Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf Inf Inf 0];
        end
    
        % Assess overfitting by comparing variance explained by 20 trials
        % relative to the rest of the dataset
        for a = 1:100
            train_select = randperm(size(Ydata,1),20);
            test_select = setdiff(1:size(Ydata,1),train_select);
            train = mean(Ydata(train_select,:),1);
            test = mean(Ydata(test_select,:),1);
            myFx = @(p) sqrt(sum((train - gammaVEP_model(xdata,p,I)).^2));
            mdl_temp = fmincon(myFx,p0,[],[],[],[],lb,ub);
            vep_fit = gammaVEP_model(xdata,mdl_temp,I);
            r = corrcoef(test,vep_fit);
            trialR2(i,I,a) = r(1,2)^2;
        end
        
        myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,I)).^2));
        mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
        [vep_fit,gamma] = gammaVEP_model(mdl_x,mdl(i,:),I);
        
        zz = 1;
        for z = 1:I
            bandwidth(i,z,:) = gamma_bandwidth(0,500,mdl(i,zz:zz+1));
            zz = zz+3;
        end
        
        [yFit(i,:)] = gammaVEP_model(xdata,mdl(i,:),I);
        r = corrcoef(ydata,yFit(i,:));
        r2(I,i) = r(1,2)^2;

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
        plot(xdata,yFit(i,:),['-' gammaC{I}],'LineWidth',2)
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
        xlabel(sprintf('r2 = %2.2f',r2(I,i)))
        title(sprintf('%2.0f',i))

        subplot(4,4,j+4)
        hold on
        for X = 1:I
             plot(mdl_x,gamma(X,:),['-' gammaC{X}])
        end
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    end

    % Determine fit on mean VEP data across individuals
    ydata = meanVEP;
    % Model guess for average model
    switch I
        case 1
            p0 = [35 75 min(ydata)];
        case 2
            p0 = [35 75 min(ydata) 25 100 max(ydata)];
        case 3
            p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata)];
        case 4
            p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata)];
        case 5
            p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata)];
        case 6
            p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata) 30 400 max(ydata)];
        case 7
            p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata) 30 400 max(ydata) 30 400 min(ydata)];
    end
    myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,I)).^2));
    Mdl = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [yFit_m,gammaM] = gammaVEP_model(xdata,Mdl,I);
    figure(205)
    plot(xdata,yFit_m,['-' gammaC{I}],'LineWidth',2)

    % plotting parameters
    figure(200)
    y = 1;
    for x = 1:(I*3)
        subplot(nGamma,3,x)
        switch x
            case [1,4,7,10,13]
                bw = bandwidth(:,y);
                hold on
                errorbar(I,mean(bw),std(bw),gammaC{I},'LineWidth',2,'MarkerFaceColor',gammaC{I},'MarkerEdgeColor',gammaC{I},'Color',gammaC{I})
                y = y+1;
            otherwise
                hold on
                errorbar(I,mean(mdl(:,x)),std(mdl(:,x)),'o','LineWidth',2,'MarkerFaceColor',gammaC{I},'MarkerEdgeColor',gammaC{I},'Color',gammaC{I})
        end
        if x == 1
            title('bandwidth')
        end
        if x == 2
            title('peak time')
        end
        if x == 3
            title('amplitude')
        end
        ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:nGamma; ax.XLim = [0 nGamma+1];
    end

    
    % Plot R2 values
    figure(202)
    hold on
    errorbar(I,mean(r2(I,:)),std(r2(I,:)),'o','LineWidth',2,'MarkerFaceColor',gammaC{I},'MarkerEdgeColor',gammaC{I},'Color',gammaC{I})
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 nGamma+1];
    
    Peak(I) = {mdl(:,2:3:end)};
    Amp(I) = {mdl(:,3:3:end)};
    Bandwidth(I) = {bandwidth};
end

% compare parameters gamma3 and gamma6
B3data = Bandwidth{1,3};
B6data = Bandwidth{1,6};
P3data = Peak{1,3};
P6data = Peak{1,6};
A3data = Amp{3};
A6data = Amp{6};

figure
subplot(3,3,1)
hold on
title('bandwidth1')
plot(B3data(:,1),B6data(:,1),'.k','MarkerSize',12)
plot([0 500],[0 500],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,4)
hold on
title('bandwidth2')
plot(B3data(:,2),B6data(:,2),'.k','MarkerSize',12)
plot([0 500],[0 500],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,7)
hold on
title('bandwidth3')
plot(B3data(:,3),B6data(:,3),'.k','MarkerSize',12)
plot([0 500],[0 500],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,2)
hold on
title('peak1')
plot(P3data(:,1),P6data(:,1),'.k','MarkerSize',12)
plot([0 500],[0 500],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,5)
hold on
title('peak2')
plot(P3data(:,2),P6data(:,2),'.k','MarkerSize',12)
plot([0 500],[0 500],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,8)
hold on
title('peak3')
plot(P3data(:,3),P6data(:,3),'.k','MarkerSize',12)
plot([0 500],[0 500],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,3)
hold on
title('amplitude1')
plot(A3data(:,1),A6data(:,1),'.k','MarkerSize',12)
plot([-2000 0],[-2000 0],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,6)
hold on
title('amplitude2')
plot(A3data(:,2),A6data(:,2),'.k','MarkerSize',12)
plot([0 4000],[0 4000],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');

subplot(3,3,9)
hold on
title('amplitude3')
plot(A3data(:,3),A6data(:,3),'.k','MarkerSize',12)
plot([-4000 0],[-4000 0],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; axis('square');