% Informed basis sets - multiple gamma functions fitting

% load train and test subjects 2021 TVST PCA paper
analysis_path = getpref('controlVEPanalysis','controlVEP_AnalysisPath');
load([analysis_path '/controlTrainTest'])

xdata = xdata.*1000; % scale time to ms

vep = [squeeze(nanmean(control_train_vep,2));squeeze(nanmean(control_test_vep,2))].*100; % x100 to correct scaling error on diopsys voltage output
subject_data = [control_train;control_test];

% truncate the part of the VEP fit to the model to improve fits
minF = 1;
maxF = 308;
vepF = vep(:,minF:maxF);
xdataF = xdata(minF:maxF);

% Determine fits on individual VEP data
yFit = zeros(size(vepF));
gamma1 = zeros(size(vepF));
gamma2 = zeros(size(vepF));
gamma3 = zeros(size(vepF));
gamma4 = zeros(size(vepF));
mdl = zeros(size(vepF,1),12);
r2 = zeros(size(vepF,1),1);
bandwidth1 = zeros(size(vepF,1),1);
bandwidth2 = zeros(size(vepF,1),1);
bandwidth3 = zeros(size(vepF,1),1);
bandwidth4 = zeros(size(vepF,1),1);

figure
j = 1;
for i = 1:size(vepF,1)
    ydata = vepF(i,:);
    myFx = @(p) sqrt(sum((ydata - gammaVEP_model4(xdataF,p)).^2));
    % Model guess
    p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata)];
    lb = [10 50 -50 10 70 1 10 90 -50 10 200 1]; 
    ub = [500 110 -1 500 150 50 500 200 -1 500 300 50];
    mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [vep_fit] = gammaVEP_model4(0:0.1:500,mdl(i,:));
    bandwidth1(i) = gamma_bandwidth(0,500,mdl(i,1:2));
    bandwidth2(i) = gamma_bandwidth(0,500,mdl(i,4:5));
    bandwidth3(i) = gamma_bandwidth(0,500,mdl(i,7:8));
    bandwidth4(i) = gamma_bandwidth(0,500,mdl(i,10:11));
    [yFit(i,:),gamma1(i,:),gamma2(i,:),gamma3(i,:),gamma4(i,:)] = gammaVEP_model4(xdataF,mdl(i,:));
    r = corrcoef(ydata,yFit(i,:));
    r2(i,:) = r(1,2)^2;
    
%     if i == 26 || i == 51
%         figure
%         j = 1;
%     else
%         j = j+1;
%     end
%     subplot(5,5,j)
    subplot(2,1,1)
    plot(xdataF,vepF(i,:),'.k')
    hold on
    plot(xdataF,yFit(i,:),'-r','LineWidth',2)
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    xlabel(sprintf('r2 = %2.2f',r2(i,:)))
    ylabel(sprintf('%2.0f',i))
    title(sprintf(['A1 = %2.0f, P1 = %2.1f, W1 = %2.0f; A2 = %2.0f, P2 = %2.1f, W2 = %2.0f;' ...
        'A3 = %2.0f, P3 = %2.1f, W3 = %2.0f; A4 = %2.0f, P4 = %2.1f, W4 = %2.0f'],...
        [mdl(i,3) mdl(i,2) bandwidth1(i) mdl(i,6) mdl(i,5) bandwidth2(i) mdl(i,9) mdl(i,8) bandwidth3(i) ...
        mdl(i,12) mdl(i,11) bandwidth4(i)]))
    hold off
    
    subplot(2,1,2)
    plot(xdataF,gamma1(i,:),'r')
    hold on
    plot(xdataF,gamma2(i,:),'b')
    plot(xdataF,gamma3(i,:),'g')
    plot(xdataF,gamma4(i,:),'m')
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    pause
    hold off
end

% Determine fit on mean VEP data across individuals
meanVEPf = mean(vepF,1);
ydata = meanVEPf;
myFx = @(p) sqrt(sum((ydata - gammaVEP_model4(xdataF,p)).^2));
% Model guess
p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata)];
Mdl = fmincon(myFx,p0,[],[],[],[],lb,ub);
[yFit_m,gamma1m,gamma2m,gamma3m,gamma4m] = gammaVEP_model4(xdataF,Mdl);
figure
plot(xdataF,meanVEPf,'.k')
hold on
plot(xdataF,yFit_m,'-r','LineWidth',2)
ax=gca; ax.TickDir = 'out'; ax.Box = 'off';

bandwidth = cat(2,bandwidth1,bandwidth2,bandwidth3,bandwidth4);

% plotting parameters
param = {'width 1','peak time 1','amplitude 1','width 2','peak time 2','amplitude 2','width 3','peak time 3','amplitude 3','width 4','peak time 4','amplitude 4'};
figure
for x = 1:12
    subplot(4,3,x)
    switch x
        case [1,4,7,10]
            tDn = mdl(:,x)./mdl(:,x-1);
            plot(ones(size(tDn)),tDn,'.k')
            hold on
            errorbar(1,mean(tDn),std(tDn),'or','LineWidth',2,'MarkerFaceColor','r')
        otherwise
            plot(ones(size(mdl(:,x))),mdl(:,x),'.k')
            hold on
            errorbar(1,mean(mdl(:,x)),std(mdl(:,x)),'or','LineWidth',2,'MarkerFaceColor','r')
    end
    ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1; ax.XTickLabel = param(x);
end



