% load train and test subjects 2021 TVST PCA paper, 6 gamma model
analysis_path = getpref('controlVEPanalysis','controlVEP_AnalysisPath');
load([analysis_path '/controlTrainTest'])

xdata = xdata.*1000; % scale time to ms

vep = [squeeze(nanmean(control_train_vep,2));squeeze(nanmean(control_test_vep,2))].*100; % x100 to correct scaling error on diopsys voltage output
meanVEP = mean(vep,1);
subject_data = [control_train;control_test];
mdl_x = 0:0.1:500;
r2 = zeros(size(vep,1),1);

gammaC = {'b','r','k','g','m','c'};
nGamma = 6;
Gamma = zeros(size(vep,1),nGamma,length(mdl_x));

gammaMdl = 'k';


figure(205)
plot(xdata,meanVEP,'.k')
hold on
ax=gca; ax.TickDir = 'out'; ax.Box = 'off';


mdl = zeros(size(vep,1),3*nGamma);
gamma = zeros(size(vep));
yFit = zeros(size(vep));
bandwidth = NaN*ones(size(vep,1),nGamma);

    % Determine fits on individual VEP data
for i = 1:size(vep,1)
    ydata = vep(i,:);
    myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,nGamma)).^2));

    p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata) 30 400 max(ydata)];
    lb = [0 0 -Inf 0 0 0 0 0 -Inf 0 0 0 0 0 -Inf 0 0 0]; 
    ub = [Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf];


    mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [vep_fit,gamma] = gammaVEP_model(mdl_x,mdl(i,:),nGamma);

    zz = 1;
    for z = 1:nGamma
        bandwidth(i,z,:) = gamma_bandwidth(0,500,mdl(i,zz:zz+1));
        zz = zz+3;
    end

    [yFit(i,:)] = gammaVEP_model(xdata,mdl(i,:),nGamma);
    r = corrcoef(ydata,yFit(i,:));
    r2(i) = r(1,2)^2;

    switch i
        case {1,9,17,25,33,41,49,57,65,73}
            figure
            j = 1;
        case {5,13,21,29,37,45,53,61,69,77}
            j = 9;
        otherwise
            j = j+1;
    end

    subplot(4,4,j)
    plot(xdata,vep(i,:),'.k')
    hold on
    plot(xdata,yFit(i,:),['-' gammaC{nGamma}],'LineWidth',2)
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    xlabel(sprintf('r2 = %2.2f',r2(i)))
    title(sprintf('subject %2.0f',i))

    subplot(4,4,j+4)
    hold on
    for X = 1:nGamma
         plot(mdl_x,gamma(X,:),['-' gammaC{X}])
    end
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    
    Gamma(i,:,:) = gamma;
end

% Determine fit on mean VEP data across individuals
ydata = meanVEP;
% Model guess for average model

p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata) 30 400 max(ydata)];
myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,nGamma)).^2));
Mdl = fmincon(myFx,p0,[],[],[],[],lb,ub);
[yFit_m,gammaM] = gammaVEP_model(xdata,Mdl,nGamma);
figure(205)
plot(xdata,yFit_m,['-' gammaC{nGamma}],'LineWidth',2)

% plotting parameters
figure(200)
subplot(1,4,1)
hold on
for x = 1:nGamma
    errorbar(x,mean(bandwidth(:,x)),std(bandwidth(:,x)),'o','LineWidth',2,'MarkerFaceColor',gammaC{x},'MarkerEdgeColor',gammaC{x},'Color',gammaC{x})
end
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:nGamma; ax.XLim = [0 nGamma+1];

subplot(1,4,4)
hold on
plot(zeros(size(r2))+rand(size(r2))-0.5,r2,'.k')
errorbar(0,mean(r2),std(r2),'o','LineWidth',2,'MarkerFaceColor',gammaMdl,'MarkerEdgeColor',gammaMdl,'Color',gammaMdl)
plot([-5 5],[1 1],'--k')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 0; ax.XLim = [-5 5]; ax.YLim = [0.7 1.1];
title('R squared')

subplot(1,4,2)
peak = mdl(:,2:3:end);
hold on
for x = 1:nGamma
    errorbar(x,mean(peak(:,x)),std(peak(:,x)),'o','LineWidth',2,'MarkerFaceColor',gammaC{x},'MarkerEdgeColor',gammaC{x},'Color',gammaC{x})
end
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:nGamma; ax.XLim = [0 nGamma+1];
title('peak time')

subplot(1,4,3)
amp = mdl(:,3:3:end);
hold on
for x = 1:nGamma
    errorbar(x,mean(amp(:,x)),std(amp(:,x)),'o','LineWidth',2,'MarkerFaceColor',gammaC{x},'MarkerEdgeColor',gammaC{x},'Color',gammaC{x})
end
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:nGamma; ax.XLim = [0 nGamma+1];
title('amplitude')

% Plot averages of gamma functions with 95% CI
figure

for x = 1:6
    subplot(1,6,x)
    hold on
    bootstat = bootstrp(1000,@mean,squeeze(Gamma(:,x,:)));
    bootstat = sort(bootstat);
    if sum(bootstat(1,:))<0
        y_dataERR1=bootstat(975,:);
        y_dataERR2=bootstat(25,:);
    else
        y_dataERR1=bootstat(25,:);
        y_dataERR2=bootstat(975,:);
    end
    
    x_ERR = cat(2,mdl_x,fliplr(mdl_x));
    y_ERR =cat(2,y_dataERR1,fliplr(y_dataERR2));
    TEMP=fill(x_ERR,y_ERR,[0.5 0.5 0.5],'EdgeColor','none');
    plot(mdl_x,bootstat(500,:),['-' gammaC{x}])
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off';
end



% Compare across age
figure
z = 1;
xdata = subject_data.age_vep;
for x = 1:6
    for y = 1:3
        subplot(6,3,z)
        hold on
        switch z
            case {1,4,7,10,13,16}
                ydata = bandwidth(:,x);
                plot(xdata,ydata,'.k','MarkerSize',12)
                title('bandwidth')
            case {2,5,8,11,14,17}
                ydata = peak(:,x);
                plot(xdata,ydata,'.k','MarkerSize',12)
                title('peak time')
            case {3,6,9,12,15,18}
                ydata = amp(:,x);
                plot(xdata,ydata,'.k','MarkerSize',12)
                title('amplitude')
        end
        lsline
        [R,P] = corrcoef(xdata,ydata);
        R2 = R(1,2)^2;
        P = P(1,2);
        ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [10 20];
        xlabel(sprintf('r2 = %2.2f, p = %2.2f',[R2 P]))
        z = z+1;
    end
end

% Linear models with age
tbl_gamma = subject_data(:,[1 448]);
tbl_gamma.bandwidth1 = bandwidth(:,1);
tbl_gamma.bandwidth2 = bandwidth(:,2);
tbl_gamma.bandwidth3 = bandwidth(:,3);
tbl_gamma.bandwidth4 = bandwidth(:,4);
tbl_gamma.bandwidth5 = bandwidth(:,5);
tbl_gamma.bandwidth6 = bandwidth(:,6);
tbl_gamma.peak1 = peak(:,1);
tbl_gamma.peak2 = peak(:,2);
tbl_gamma.peak3 = peak(:,3);
tbl_gamma.peak4 = peak(:,4);
tbl_gamma.peak5 = peak(:,5);
tbl_gamma.peak6 = peak(:,6);
tbl_gamma.amp1 = amp(:,1);
tbl_gamma.amp2 = amp(:,2);
tbl_gamma.amp3 = amp(:,3);
tbl_gamma.amp4 = amp(:,4);
tbl_gamma.amp5 = amp(:,5);
tbl_gamma.amp6 = amp(:,6);

mdl_gamma = fitglm(tbl_gamma,...
    'age_vep ~ bandwidth1 + peak1 + amp1 + bandwidth2 + peak2 + amp2 + bandwidth3 + peak3 + amp3 + bandwidth4 + peak4 + amp4 + bandwidth5 + peak5 + amp5 + bandwidth6 + peak6 + amp6');

mdl_gamma2 = fitglm(tbl_gamma,'age_vep ~ bandwidth1 + peak1 + amp1 + bandwidth2 + peak2 + amp2');
