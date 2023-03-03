% Informed basis sets - multiple gamma functions fitting

% VEP to test model on

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');

load([data_path '/cleaned_VEP.mat'])

% Get the mean VEP per participant
xdata = cell2mat(control_vep(1,3)).*1000;

vep = zeros(length(control_vep),length(xdata));
for i = 1:length(control_vep)
    temp = cell2mat(control_vep(i,2));
    vep (i,:) = mean(temp,1).*10;
end

% bad_fits = [44 45 46 47 48 49 51 53 54 60];
% vep = vep(bad_fits,:);

% truncate the part of the VEP fit to the model to improve fits
minF = 1;
maxF = 308;
vepF = vep(:,minF:maxF);
xdataF = xdata(minF:maxF);

% Determine fits on individual VEP data
yFit = zeros(size(vepF));
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
    myFx = @(p) sqrt(sum((ydata - gammaVEP_model3(xdataF,p)).^2));
    % Model guess
    p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata)];
    lb = [10 50 -50 10 70 0 10 90 -50 10 200 0]; 
    ub = [500 110 0 500 150 50 500 200 0 500 300 50];
    mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [vep_fit] = gammaVEP_model3(0:0.1:500,mdl(i,:));
    bandwidth1(i) = gamma_bandwidth(0,500,mdl(i,1:2));
    bandwidth2(i) = gamma_bandwidth(0,500,mdl(i,4:5));
    bandwidth3(i) = gamma_bandwidth(0,500,mdl(i,7:8));
    bandwidth4(i) = gamma_bandwidth(0,500,mdl(i,10:11));
    yFit(i,:) = gammaVEP_model3(xdataF,mdl(i,:));
    r = corrcoef(ydata,yFit(i,:));
    r2(i,:) = r(1,2)^2;
    
%     if i == 26 || i == 51
%         figure
%         j = 1;
%     else
%         j = j+1;
%     end
%     subplot(5,5,j)
    
    plot(xdata,vep(i,:),'.k')
    hold on
    plot(xdataF,yFit(i,:),'-r','LineWidth',2)
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    xlabel(sprintf('r2 = %2.2f',r2(i,:)))
    ylabel(sprintf('%2.0f',i))
    title(sprintf(['A1 = %2.0f, P1 = %2.1f, W1 = %2.0f; A2 = %2.0f, P2 = %2.1f, W2 = %2.0f;' ...
        'A3 = %2.0f, P3 = %2.1f, W3 = %2.0f; A4 = %2.0f, P4 = %2.1f, W4 = %2.0f'],...
        [mdl(i,3) mdl(i,2) bandwidth1(i) mdl(i,6) mdl(i,5) bandwidth2(i) mdl(i,9) mdl(i,8) bandwidth3(i) ...
        mdl(i,12) mdl(i,11) bandwidth4(i)]))
    pause
    hold off
end

% Determine fit on mean VEP data across individuals
meanVEP = mean(vep,1);
meanVEPf = meanVEP(:,minF:maxF);
ydata = meanVEPf;
myFx = @(p) sqrt(sum((ydata - gammaVEP_model3(xdataF,p)).^2));
% Model guess
p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata)];
Mdl = fmincon(myFx,p0,[],[],[],[],lb,ub);
yFit_m = gammaVEP_model3(xdataF,Mdl);
figure
plot(xdata,meanVEP,'.k')
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



