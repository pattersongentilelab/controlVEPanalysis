% Use single gamma model (can adjust to different number of gammas)


data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])
saveFig_path = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/';

%% select control subjects with two sessions within 6 months of each other
control_subject = find(cleaned_vep_files.subjecttype=='Control');
cleaned_vep = cleaned_vep(control_subject,:);
control_vep_subjects = cleaned_vep_files(control_subject,:);
xdata = cleaned_vep{1,3}.*1000; %convert to ms

unique_ID = unique(control_vep_subjects.uniqueID);
vep = zeros(90,2,length(xdata));
subject_sess1 = control_vep_subjects(1,:);
subject_sess2 = control_vep_subjects(1,:);

counter=1;
counter2=1;
for x=1:length(unique_ID)
    temp_loc = find(cell2mat(cleaned_vep(:,1))==unique_ID(x,:));
    if length(temp_loc)>1 % select only subjects with multiple sessions <6 months apart
        age_range = control_vep_subjects.age_vep(temp_loc(1:2),:);
        if diff(age_range,[],1) < 0.5
        temp_loc = temp_loc(1:2); % select only the first two sessions
        temp_ydata = cleaned_vep(temp_loc,4);
            for y=1:size(temp_ydata,1)
                vep(counter2,y,:) = nanmean(cell2mat(temp_ydata(y,:))).*100; %correct amplification issue with diopsys
                counter=counter+1;
            end
        subject_sess1(counter2,:) = control_vep_subjects(temp_loc(1),:);
        subject_sess2(counter2,:) = control_vep_subjects(temp_loc(2),:);
        counter2 = counter2+1;
        end
    end
    clear temp_loc temp_ydata
end

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
nGamma = 4;

xdata = xdata(:,1:xdata_end);
vep = vep(:,:,1:xdata_end);

meanVEP = squeeze(mean(vep,1));

r_val = zeros(size(vep,1),2);
r2 = zeros(size(vep,1),2);

%% Check for local minima and maxima for model guesses

Amp75 = zeros(size(vep,1),2);
Peak75 = zeros(size(vep,1),2);
Bw75 = zeros(size(vep,1),2);
Amp100 = zeros(size(vep,1),2);
Peak100 = zeros(size(vep,1),2);
Bw100 = zeros(size(vep,1),2);
Amp135 = zeros(size(vep,1),2);
Peak135 = zeros(size(vep,1),2);
Bw135 = zeros(size(vep,1),2);
Amp220 = zeros(size(vep,1),2);
Peak220 = zeros(size(vep,1),2);
Bw220 = zeros(size(vep,1),2);

for i = 1:size(vep,1)
    for j = 1:2
        ydata = squeeze(vep(i,j,:));
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

        Amp75(i,j) = amp75;
        Peak75(i,j) = peak75;
        Bw75(i,j) = bw75;
        Amp100(i,j) = amp100;
        Peak100(i,j) = peak100;
        Bw100(i,j) = bw100;
        Amp135(i,j) = amp135;
        Peak135(i,j) = peak135;
        Bw135(i,j) = bw135;
        Amp220(i,j) = amp220;
        Peak220(i,j) = peak220;
        Bw220(i,j) = bw220;
    end
end

%% Determine fits on individual VEP data
mdl = zeros(size(vep,1),2,3*nGamma);
Gamma = zeros(size(vep,1),2,nGamma,length(mdl_x));
yFit = zeros(size(vep));
bandwidth = NaN*ones(size(vep,1),2,nGamma);
sess_marker = {'o','s'};

for I = 1:2 % loop across session

    for i = 1:size(vep,1)

        ydata = squeeze(vep(i,I,:))'; % averaged across trials, already corrected for diopsys amplification error above
        p0 = [Bw75(i,j) Peak75(i,j) Amp75(i,j) Bw100(i,j) Peak100(i,j) Amp100(i,j) Bw135(i,j) Peak135(i,j) Amp135(i,j) Bw220(i,j) Peak220(i,j) Amp220(i,j)];
        lb = [max([30 Bw75(i,j)-5]) Peak75(i,j)-2 Amp75(i,j)*1.1 max([20 Bw100(i,j)-5]) Peak100(i,j)-3 Amp100(i,j)*0.9 max([15 Bw135(i,j)-5]) Peak135(i,j)-5 Amp135(i,j)*1.1 max([15 Bw220(i,j)-5]) Peak220(i,j)-5 Amp220(i,j)*0.9]; 
        ub = [min([110 Bw75(i,j)+5]) Peak75(i,j)+2 Amp75(i,j)*0.9 min([110 Bw100(i,j)+5]) Peak100(i,j)+3 Amp100(i,j)*1.1 min([110 Bw135(i,j)+5]) Peak135(i,j)+5 -Amp135(i,j)*0.9 min([100 Bw220(i,j)+5]) Peak220(i,j)+5 Amp220(i,j)*1.1];
       

        myFx = @(p) sqrt(sum((ydata - gammaVEP_model(xdata,p,nGamma)).^2));
        mdl(i,I,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
        [vep_fit,gamma] = gammaVEP_model(mdl_x,squeeze(mdl(i,I,:)),nGamma);
        Gamma(i,I,:,:) = gamma;

        for z = 1:nGamma
            bandwidth(i,I,z) = gamma_bandwidth(mdl_x,gamma(z,:));
        end

        [yFit(i,I,:)] = gammaVEP_model(xdata,squeeze(mdl(i,I,:)),nGamma);
        r = corrcoef(ydata,squeeze(yFit(i,I,:)));
        r_val(i,I) = r(1,2);
        r2(i,I) = r(1,2)^2;
    end
    % Plot R values
    figure(202)
    hold on
    plot(ones(size(r_val(:,I))),r_val(:,I),'.','Color',[0.5 0.5 0.5])
    errorbar(1,mean(r_val(:,I)),std(r_val(:,I)),sess_marker{I},'LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [0 2];
end

peak = mdl(:,:,2:3:end);
amp = mdl(:,:,3:3:end);

%% Compare session 1 and session 2 parameters

% Informed basis set
fig = figure(300);
z = 1;
for x = 1:3
    for y = 1:nGamma
        subplot(3,nGamma,z)
        hold on
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    
        switch x
            case 1
                plot(squeeze(bandwidth(:,1,y)),squeeze(bandwidth(:,2,y)),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})
                plot([0 120],[0 120],'--')
                [rr,pp] = corrcoef(squeeze(bandwidth(:,1,y)),squeeze(bandwidth(:,2,y)));
                title(sprintf('bandwidth %1d',y))
                ylim([0 200])
                xlim([0 200])
            case 2
                plot(squeeze(peak(:,1,y)),squeeze(peak(:,2,y)),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})
                plot([0 400],[0 400],'--')
                plot([10 20],[lb((y*3)-1) lb((y*3)-1)],'--')
                plot([10 20],[ub((y*3)-1) ub((y*3)-1)],'--')
                [rr,pp] = corrcoef(squeeze(peak(:,1,y)),squeeze(peak(:,2,y)));
                title(sprintf('peak time %1d',y))
                ylim([0 400])
                xlim([0 400])
            case 3
                plot(squeeze(abs(amp(:,1,y))),squeeze(abs(amp(:,2,y))),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor','w','Color',gammaC{y})
                plot([0 max(max(max(amp)))],[0 max(max(max(amp)))],'--')
                [rr,pp] = corrcoef(squeeze(abs(amp(:,1,y))),squeeze(abs(amp(:,2,y))));
                title(sprintf('amplitude %1d',y))
                ylim([0 max(max(max(amp)))])
                xlim([0 max(max(max(amp)))])
        end
        xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))
        z = z+1;
    end
end

fig_name = [saveFig_path 'twoSessionParamsCorr'];
print(fig,fig_name,'-dpdf','-painters')

%% Bland-Altman plot

% Informed basis set

fig = figure(302);
subplot(3,4,1)
diff_n75bw_sess = diff(squeeze(bandwidth(:,:,1)),[],2);
mean_n75bw_sess = mean(squeeze(bandwidth(:,:,1)),2);
hold on
plot(mean_n75bw_sess,diff_n75bw_sess,'o','MarkerFaceColor','b','MarkerEdgeColor','w','Color','b')
plot([min(mean_n75bw_sess) max(mean_n75bw_sess)],[mean(diff_n75bw_sess) mean(diff_n75bw_sess)],'--k')
plot([min(mean_n75bw_sess) max(mean_n75bw_sess)],[mean(diff_n75bw_sess)+std(diff_n75bw_sess)*1.96 mean(diff_n75bw_sess)+std(diff_n75bw_sess)*1.96],'--k')
plot([min(mean_n75bw_sess) max(mean_n75bw_sess)],[mean(diff_n75bw_sess)+std(diff_n75bw_sess)*-1.96 mean(diff_n75bw_sess)+std(diff_n75bw_sess)*-1.96],'--k')
text(max(mean_n75bw_sess)-2,mean(diff_n75bw_sess)+std(diff_n75bw_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n75bw_sess)+std(diff_n75bw_sess)*-1.96));
text(max(mean_n75bw_sess)-2,mean(diff_n75bw_sess)+std(diff_n75bw_sess)*1.96+2,sprintf('%2.1f',mean(diff_n75bw_sess)+std(diff_n75bw_sess)*-1.96));
title('Bland Altman N75 bandwidth')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75bw_sess), max(mean_n75bw_sess)]; ax.YLim = [-20,20];

subplot(3,4,2)
diff_p100bw_sess = diff(squeeze(bandwidth(:,:,2)),[],2);
mean_p100bw_sess = mean(squeeze(bandwidth(:,:,2)),2);
hold on
plot(mean_p100bw_sess,diff_p100bw_sess,'o','MarkerFaceColor','r','MarkerEdgeColor','w','Color','r')
plot([min(mean_p100bw_sess) max(mean_p100bw_sess)],[mean(diff_p100bw_sess) mean(diff_p100bw_sess)],'--k')
plot([min(mean_p100bw_sess) max(mean_p100bw_sess)],[mean(diff_p100bw_sess)+std(diff_p100bw_sess)*1.96 mean(diff_p100bw_sess)+std(diff_p100bw_sess)*1.96],'--k')
plot([min(mean_p100bw_sess) max(mean_p100bw_sess)],[mean(diff_p100bw_sess)+std(diff_p100bw_sess)*-1.96 mean(diff_p100bw_sess)+std(diff_p100bw_sess)*-1.96],'--k')
text(max(mean_p100bw_sess)-2,mean(diff_p100bw_sess)+std(diff_p100bw_sess)*-1.96-2,sprintf('%2.1f',mean(diff_p100bw_sess)+std(diff_p100bw_sess)*-1.96));
text(max(mean_p100bw_sess)-2,mean(diff_p100bw_sess)+std(diff_p100bw_sess)*1.96+2,sprintf('%2.1f',mean(diff_p100bw_sess)+std(diff_p100bw_sess)*1.96));

title('Bland Altman P100 bandwidth')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_p100bw_sess), max(mean_p100bw_sess)]; ax.YLim = [-20,20];

subplot(3,4,3)
diff_n135bw_sess = diff(squeeze(bandwidth(:,:,3)),[],2);
mean_n135bw_sess = mean(squeeze(bandwidth(:,:,3)),2);
hold on
plot(mean_n135bw_sess,diff_n135bw_sess,'o','MarkerFaceColor','m','MarkerEdgeColor','w','Color','m')
plot([min(mean_n135bw_sess) max(mean_n135bw_sess)],[mean(diff_n135bw_sess) mean(diff_n135bw_sess)],'--k')
plot([min(mean_n135bw_sess) max(mean_n135bw_sess)],[mean(diff_n135bw_sess)+std(diff_n135bw_sess)*1.96 mean(diff_n135bw_sess)+std(diff_n135bw_sess)*1.96],'--k')
plot([min(mean_n135bw_sess) max(mean_n135bw_sess)],[mean(diff_n135bw_sess)+std(diff_n135bw_sess)*-1.96 mean(diff_n135bw_sess)+std(diff_n135bw_sess)*-1.96],'--k')
text(max(mean_n135bw_sess)-2,mean(diff_n135bw_sess)+std(diff_n135bw_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n135bw_sess)+std(diff_n135bw_sess)*-1.96));
text(max(mean_n135bw_sess)-2,mean(diff_n135bw_sess)+std(diff_n135bw_sess)*1.96+2,sprintf('%2.1f',mean(diff_n135bw_sess)+std(diff_n135bw_sess)*1.96));
title('Bland Altman N135 bandwidth')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n135bw_sess), max(mean_n135bw_sess)]; ax.YLim = [-20,20];

subplot(3,4,4)
diff_pLatebw_sess = diff(squeeze(bandwidth(:,:,4)),[],2);
mean_pLatebw_sess = mean(squeeze(bandwidth(:,:,4)),2);
hold on
plot(mean_pLatebw_sess,diff_pLatebw_sess,'o','MarkerFaceColor','g','MarkerEdgeColor','w','Color','g')
plot([min(mean_pLatebw_sess) max(mean_pLatebw_sess)],[mean(diff_pLatebw_sess) mean(diff_pLatebw_sess)],'--k')
plot([min(mean_pLatebw_sess) max(mean_pLatebw_sess)],[mean(diff_pLatebw_sess)+std(diff_pLatebw_sess)*1.96 mean(diff_pLatebw_sess)+std(diff_pLatebw_sess)*1.96],'--k')
plot([min(mean_pLatebw_sess) max(mean_pLatebw_sess)],[std(diff_pLatebw_sess)*-1.96 std(diff_pLatebw_sess)*-1.96],'--k')
text(max(mean_pLatebw_sess)-2,mean(diff_pLatebw_sess)+std(diff_pLatebw_sess)*-1.96-2,sprintf('%2.1f',mean(diff_pLatebw_sess)+std(diff_pLatebw_sess)*-1.96));
text(max(mean_pLatebw_sess)-2,mean(diff_pLatebw_sess)+std(diff_pLatebw_sess)*1.96+2,sprintf('%2.1f',mean(diff_pLatebw_sess)+std(diff_pLatebw_sess)*1.96));
title('Bland Altman Late bandwidth')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_pLatebw_sess), max(mean_pLatebw_sess)]; ax.YLim = [-20,20];


subplot(3,4,5)
diff_n75pk_sess = diff(squeeze(peak(:,:,1)),[],2);
mean_n75pk_sess = mean(squeeze(peak(:,:,1)),2);
hold on
plot(mean_n75pk_sess,diff_n75pk_sess,'o','MarkerFaceColor','b','MarkerEdgeColor','w','Color','b')
plot([min(mean_n75pk_sess) max(mean_n75pk_sess)],[mean(diff_n75pk_sess) mean(diff_n75pk_sess)],'--k')
plot([min(mean_n75pk_sess) max(mean_n75pk_sess)],[mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96 mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96],'--k')
plot([min(mean_n75pk_sess) max(mean_n75pk_sess)],[mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96 mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96],'--k')
text(max(mean_n75pk_sess)-2,mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n75pk_sess)+std(diff_n75pk_sess)*-1.96));
text(max(mean_n75pk_sess)-2,mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96+2,sprintf('%2.1f',mean(diff_n75pk_sess)+std(diff_n75pk_sess)*1.96));
title('Bland Altman N75 peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75pk_sess), max(mean_n75pk_sess)]; ax.YLim = [-12,12];

subplot(3,4,6)
diff_p100pk_sess = diff(squeeze(peak(:,:,2)),[],2);
mean_p100pk_sess = mean(squeeze(peak(:,:,2)),2);
hold on
plot(mean_p100pk_sess,diff_p100pk_sess,'o','MarkerFaceColor','r','MarkerEdgeColor','w','Color','r')
plot([min(mean_p100pk_sess) max(mean_p100pk_sess)],[mean(diff_p100pk_sess) mean(diff_p100pk_sess)],'--k')
plot([min(mean_p100pk_sess) max(mean_p100pk_sess)],[mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96 mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96],'--k')
plot([min(mean_p100pk_sess) max(mean_p100pk_sess)],[mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96 mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96],'--k')
text(max(mean_p100pk_sess)-2,mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_p100pk_sess)+std(diff_p100pk_sess)*-1.96));
text(max(mean_p100pk_sess)-2,mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96+2,sprintf('%2.1f',mean(diff_p100pk_sess)+std(diff_p100pk_sess)*1.96));
title('Bland Altman P100 peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_p100pk_sess), max(mean_p100pk_sess)]; ax.YLim = [-12,12];

subplot(3,4,7)
diff_n135pk_sess = diff(squeeze(peak(:,:,3)),[],2);
mean_n135pk_sess = mean(squeeze(peak(:,:,3)),2);
hold on
plot(mean_n135pk_sess,diff_n135pk_sess,'o','MarkerFaceColor','m','MarkerEdgeColor','w','Color','m')
plot([min(mean_n135pk_sess) max(mean_n135pk_sess)],[mean(diff_n135pk_sess) mean(diff_n135pk_sess)],'--k')
plot([min(mean_n135pk_sess) max(mean_n135pk_sess)],[mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96 mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96],'--k')
plot([min(mean_n135pk_sess) max(mean_n135pk_sess)],[mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96 mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96],'--k')
text(max(mean_n135pk_sess)-2,mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n135pk_sess)+std(diff_n135pk_sess)*-1.96));
text(max(mean_n135pk_sess)-2,mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96+2,sprintf('%2.1f',mean(diff_n135pk_sess)+std(diff_n135pk_sess)*1.96));
title('Bland Altman N135 peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n135pk_sess), max(mean_n135pk_sess)]; ax.YLim = [-12,12];

subplot(3,4,8)
diff_pLatepk_sess = diff(squeeze(peak(:,:,4)),[],2);
mean_pLatepk_sess = mean(squeeze(peak(:,:,4)),2);
hold on
plot(mean_pLatepk_sess,diff_pLatepk_sess,'o','MarkerFaceColor','g','MarkerEdgeColor','w','Color','g')
plot([min(mean_pLatepk_sess) max(mean_pLatepk_sess)],[mean(diff_pLatepk_sess) mean(diff_pLatepk_sess)],'--k')
plot([min(mean_pLatepk_sess) max(mean_pLatepk_sess)],[mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96 mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96],'--k')
plot([min(mean_pLatepk_sess) max(mean_pLatepk_sess)],[mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96 mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96],'--k')
text(max(mean_pLatepk_sess)-2,mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96-2,sprintf('%2.1f',mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*-1.96));
text(max(mean_pLatepk_sess)-2,mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96+2,sprintf('%2.1f',mean(diff_pLatepk_sess)+std(diff_pLatepk_sess)*1.96));
title('Bland Altman Late peak')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_pLatepk_sess), max(mean_pLatepk_sess)]; ax.YLim = [-12,12];


subplot(3,4,9)
diff_n75am_sess = diff(squeeze(amp(:,:,1)),[],2);
mean_n75am_sess = mean(squeeze(amp(:,:,1)),2);
hold on
plot(mean_n75am_sess,diff_n75am_sess,'o','MarkerFaceColor','b','MarkerEdgeColor','w','Color','b')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess) mean(diff_n75am_sess)],'--k')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96 mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96],'--k')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96 mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96],'--k')
text(max(mean_n75am_sess)-2,mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96));
text(max(mean_n75am_sess)-2,mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96+2,sprintf('%2.1f',mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96));
title('Bland Altman N75 amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75am_sess), max(mean_n75am_sess)]; ax.YLim = [-6,6];

subplot(3,4,10)
diff_p100am_sess = diff(squeeze(amp(:,:,2)),[],2);
mean_p100am_sess = mean(squeeze(amp(:,:,2)),2);
hold on
plot(mean_p100am_sess,diff_p100am_sess,'o','MarkerFaceColor','r','MarkerEdgeColor','w','Color','r')
plot([min(mean_p100am_sess) max(mean_p100am_sess)],[mean(diff_p100am_sess) mean(diff_p100am_sess)],'--k')
plot([min(mean_p100am_sess) max(mean_p100am_sess)],[mean(diff_p100am_sess)+std(diff_p100am_sess)*1.96 mean(diff_p100am_sess)+std(diff_p100am_sess)*1.96],'--k')
plot([min(mean_p100am_sess) max(mean_p100am_sess)],[mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96 mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96],'--k')
text(max(mean_p100am_sess)-2,mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96-2,sprintf('%2.1f',mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96));
text(max(mean_p100am_sess)-2,mean(diff_p100am_sess)+std(diff_p100am_sess)*-1.96+2,sprintf('%2.1f',mean(diff_p100am_sess)+std(diff_p100am_sess)*1.96));
title('Bland Altman P100 amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_p100am_sess), max(mean_p100am_sess)]; ax.YLim = [-6,6];

subplot(3,4,11)
diff_n135am_sess = diff(squeeze(amp(:,:,3)),[],2);
mean_n135am_sess = mean(squeeze(amp(:,:,3)),2);
hold on
plot(mean_n135am_sess,diff_n135am_sess,'o','MarkerFaceColor','m','MarkerEdgeColor','w','Color','m')
plot([min(mean_n135am_sess) max(mean_n135am_sess)],[mean(diff_n135am_sess) mean(diff_n135am_sess)],'--k')
plot([min(mean_n135am_sess) max(mean_n135am_sess)],[mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96 mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96],'--k')
plot([min(mean_n135am_sess) max(mean_n135am_sess)],[mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96 mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96],'--k')
text(max(mean_n135am_sess)-2,mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96-2,sprintf('%2.1f',mean(diff_n135am_sess)+std(diff_n135am_sess)*-1.96));
text(max(mean_n135am_sess)-2,mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96+2,sprintf('%2.1f',mean(diff_n135am_sess)+std(diff_n135am_sess)*1.96));
title('Bland Altman N135 amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n135am_sess), max(mean_n135am_sess)]; ax.YLim = [-6,6];

subplot(3,4,12)
diff_pLateam_sess = diff(squeeze(amp(:,:,4)),[],2);
mean_pLateam_sess = mean(squeeze(amp(:,:,4)),2);
hold on
plot(mean_pLateam_sess,diff_pLateam_sess,'o','MarkerFaceColor','g','MarkerEdgeColor','w','Color','g')
plot([min(mean_pLateam_sess) max(mean_pLateam_sess)],[mean(diff_pLateam_sess) mean(diff_pLateam_sess)],'--k')
plot([min(mean_pLateam_sess) max(mean_pLateam_sess)],[mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96 mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96],'--k')
plot([min(mean_pLateam_sess) max(mean_pLateam_sess)],[mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96 mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96],'--k')
text(max(mean_pLateam_sess)-2,mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96-2,sprintf('%2.1f',mean(diff_pLateam_sess)+std(diff_pLateam_sess)*-1.96));
text(max(mean_pLateam_sess)-2,mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96+2,sprintf('%2.1f',mean(diff_pLateam_sess)+std(diff_pLateam_sess)*1.96));
title('Bland Altman Late amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_pLateam_sess), max(mean_pLateam_sess)]; ax.YLim = [-6,6];

fig_name = [saveFig_path 'twoSessionParamsBlandAltman'];
print(fig,fig_name,'-dpdf','-painters')

%% Direct comparison Bland Altman plots to standard ICSEV analysis

subplot(3,4,9)
diff_n75am_sess = diff(squeeze(amp(:,:,1)),[],2);
mean_n75am_sess = mean(squeeze(amp(:,:,1)),2);
hold on
plot(mean_n75am_sess,diff_n75am_sess,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5])
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess) mean(diff_n75am_sess)],'--k')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96 mean(diff_n75am_sess)+std(diff_n75am_sess)*1.96],'--k')
plot([min(mean_n75am_sess) max(mean_n75am_sess)],[mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96 mean(diff_n75am_sess)+std(diff_n75am_sess)*-1.96],'--k')
title('Bland Altman N75 amplitude')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75am_sess), max(mean_n75am_sess)]; ax.YLim = [diff([max(mean_n75am_sess) min(mean_n75am_sess)]),diff([min(mean_n75am_sess) max(mean_n75am_sess)])];


%% Plot individual VEP to check fits

figure
for i = 1:size(vep,1)
    for j = 1:2
        subplot(1,2,j)
        hold on
        plot(mdl_x,squeeze(squeeze(sum(Gamma(i,j,:,:),3))),'c')
        plot(xdata,squeeze(vep(i,j,:)),'Color',[0.5 0.5 0.5])
        for X = 1:nGamma
             plot(mdl_x,squeeze(squeeze(Gamma(i,j,X,:))),['-' gammaC{X}])
             text(squeeze(peak(i,j,1)),squeeze(amp(i,j,1)),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[squeeze(bandwidth(i,j,1)) squeeze(peak(i,j,1)) squeeze(amp(i,j,1))]));
             text(squeeze(peak(i,j,2)),squeeze(amp(i,j,2)),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[squeeze(bandwidth(i,j,2)) squeeze(peak(i,j,2)) squeeze(amp(i,j,2))]));
             text(squeeze(peak(i,j,3)),squeeze(amp(i,j,3)),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[squeeze(bandwidth(i,j,3)) squeeze(peak(i,j,3)) squeeze(amp(i,j,3))]));
             text(squeeze(peak(i,j,4)),squeeze(amp(i,j,4)),sprintf('bw = %2.2f \n pt = %2.2f \n amp = %2.2f',[squeeze(bandwidth(i,j,4)) squeeze(peak(i,j,4)) squeeze(amp(i,j,4))]));        
        end
        
        plot(Peak75(i,j),Amp75(i,j),'+b')
        text(Peak75(i,j)+0.1,Amp75(i,j)+0.1,sprintf('bw guess = %2.2f',Bw75(i,j)));
        plot(Peak100(i,j),Amp100(i,j),'+r')
        text(Peak100(i,j)+0.1,Amp100(i,j)+0.1,sprintf('bw guess = %2.2f',Bw100(i,j)));
        plot(Peak135(i,j),Amp135(i,j),'+m')
        text(Peak135(i,j)+0.1,Amp135(i,j)+0.1,sprintf('bw guess = %2.2f',Bw135(i,j)));
        plot(Peak220(i,j),Amp220(i,j),'+g')
        text(Peak220(i,j)+0.1,Amp220(i,j)+0.1,sprintf('bw guess = %2.2f',Bw220(i,j)));
        title(sprintf('diff bw n75 = %2.2f, p100 = %2.2f, n135 = %2.2f, late = %2.2f',[diff_n75bw_sess(i) diff_p100bw_sess(i) diff_n135bw_sess(i) diff_pLatebw_sess(i)]));   
        xlabel(sprintf('subject %2.0f, Fit r = %2.2f',[i r_val(i,j)]))
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.YLim = [-45 45]; ax.XLim = [0 time_end];
    end
    pause(1)
    subplot(1,2,1)
    clf
    subplot(1,2,2)
    clf
end
