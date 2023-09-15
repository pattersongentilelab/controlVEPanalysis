% Use gamma model with 4 components, more fixed variables to account for the strong 
% correlation between peak time and preceding bandwidth


data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])
save_indivVEP_path = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/indivVEP_wFits/';
load([data_path '/neuroActive_meds.mat']) 

%% Select control subjects, first session, who do not have a history of migraine (med hx1), chronic headache (med hx2), vision or reading therapy (med hx11+12) or strabismus/amblyopia (15 - 18)
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

% %remove participants on neuroactive meds
% temp = setdiff(subject.uniqueID,neuro_active_meds);
% subject = subject(ismember(subject.uniqueID,temp),:);
% vep = vep(ismember(subject.uniqueID,temp),:);


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

r_val = zeros(size(vep,1),1);
r_val300 = zeros(size(vep,1),1);


%% Determine fits on individual VEP data
mdl = zeros(size(vep,1),9);
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
    
    p0 = [Bw75(i,:) Peak75(i,:) Amp75(i,:) Bw100(i,:) Amp100(i,:) Bw135(i,:) Amp135(i,:) Bw220(i,:) Amp220(i,:)];
    lb = [max([30 Bw75(i,:)-5]) Peak75(i,:)-2 Amp75(i,:)*1.1 max([20 Bw100(i,:)-5]) Amp100(i,:)*0.9 max([15 Bw135(i,:)-5]) Amp135(i,:)*1.1 max([15 Bw220(i,:)-5]) Amp220(i,:)*0.9]; 
    ub = [min([110 Bw75(i,:)+5]) Peak75(i,:)+2 Amp75(i,:)*0.9 min([110 Bw100(i,:)+5]) Amp100(i,:)*1.1 min([110 Bw135(i,:)+5]) -Amp135(i,:)*0.9 min([100 Bw220(i,:)+5]) Amp220(i,:)*1.1];

    myFx = @(p) sqrt(sum((ydata - gammaVEP_modelReduced(xdata,p)).^2));
    mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [vep_fit,gamma] = gammaVEP_modelReduced(mdl_x,mdl(i,:));
    Gamma(i,:,:) = gamma;
    
    for z = 1:nGamma
        bandwidth(i,z,:) = gamma_bandwidth(mdl_x,gamma(z,:));
    end

    [yFit(i,:)] = gammaVEP_modelReduced(xdata,mdl(i,:));
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

    
% Plot R values
figure(202)
hold on
plot(ones(size(r_val)),r_val,'.','Color',[0.5 0.5 0.5])
errorbar(1,mean(r_val),std(r_val),'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
plot(2*ones(size(r_val300)),r_val300,'.','Color',[0.5 0.5 0.5])
errorbar(2,mean(r_val300),std(r_val300),'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 3];




