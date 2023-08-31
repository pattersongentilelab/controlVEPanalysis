% Use single gamma model (can adjust to different number of gammas)


data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');
load([data_path '/randControlTrainTest.mat']) % participants selected for train and test
load([data_path '/cleaned_VEP.mat'])
save_indivVEP_path = '/Users/pattersonc/Library/CloudStorage/OneDrive-Children''sHospitalofPhiladelphia/Research/Minds Matter/Figures/InformedBasisSet/indivVEP_wFits/';

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

%% compile peak analysis
lat_n75 = [subject_sess1.Left_Cursor_Lat subject_sess2.Left_Cursor_Lat];
lat_p100 = [subject_sess1.Right_Cursor_Lat subject_sess2.Right_Cursor_Lat];

peak_n75 = [subject_sess1.Left_Cursor_Amp subject_sess2.Left_Cursor_Amp];
peak_p100 = [subject_sess1.Right_Cursor_Amp subject_sess2.Right_Cursor_Amp];

peak_diff1 = abs(diff([peak_n75(:,1) peak_p100(:,1)],[],2));
peak_diff2 = abs(diff([peak_n75(:,2) peak_p100(:,2)],[],2));
peak_diff = [peak_diff1 peak_diff2];

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

%% Determine fits on individual VEP data
mdl = zeros(size(vep,1),2,3*nGamma);
Gamma = zeros(size(vep,1),2,nGamma,length(mdl_x));
yFit = zeros(size(vep));
bandwidth = NaN*ones(size(vep,1),2,nGamma);
sess_marker = {'o','s'};

for I = 1:2 % loop across session

    for i = 1:size(vep,1)

        ydata = squeeze(vep(i,I,:))'; % averaged across trials, already corrected for diopsys amplification error above

        p0 = [35 75 min(ydata) 40 100 max(ydata) 40 135 min(ydata) 25 250 max(ydata)];
        lb = [30 60 min(ydata)*1.2 30 90 0.5 30 120 min(ydata)*1.2 10 180 0.5]; 
        ub = [60 90 -0.5 150 120 max(ydata)*1.2 150 170 -0.5 50 300 max(ydata)*1.2];

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
figure(300)
z = 1;
for x = 1:3
    for y = 1:nGamma
        subplot(3,nGamma,z)
        hold on
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    
        switch x
            case 1
                plot(squeeze(bandwidth(:,1,y)),squeeze(bandwidth(:,2,y)),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
                plot([0 120],[0 120],'--')
                [rr,pp] = corrcoef(squeeze(bandwidth(:,1,y)),squeeze(bandwidth(:,2,y)));
                title(sprintf('bandwidth %1d',y))
                ylim([0 200])
            case 2
                plot(squeeze(peak(:,1,y)),squeeze(peak(:,2,y)),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
                plot([0 400],[0 400],'--')
                plot([10 20],[lb((y*3)-1) lb((y*3)-1)],'--')
                plot([10 20],[ub((y*3)-1) ub((y*3)-1)],'--')
                [rr,pp] = corrcoef(squeeze(peak(:,1,y)),squeeze(peak(:,2,y)));
                title(sprintf('peak time %1d',y))
                ylim([0 400])
                xlim([0 400])
            case 3
                plot(squeeze(abs(amp(:,1,y))),squeeze(abs(amp(:,2,y))),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
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


% peak analysis

figure(301)
subplot(1,3,1)
hold on
plot(peak_diff(:,1),peak_diff(:,2),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
plot([0 max(max(peak_diff))],[0 max(max(peak_diff))],'--')
[rr,pp] = corrcoef(peak_diff(:,1),peak_diff(:,2));
title('N75-P100 peak difference')
ylim([0 max(max(peak_diff))])
xlim([0 max(max(peak_diff))])
xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))

subplot(1,3,2)
hold on
plot(lat_n75(:,1),lat_n75(:,2),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
plot([0 400],[0 400],'--')
[rr,pp] = corrcoef(lat_n75(:,1),lat_n75(:,2));
title('N75 peak time')
ylim([0 400])
xlim([0 400])
xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))

subplot(1,3,3)
hold on
plot(lat_p100(:,1),lat_p100(:,2),'o','MarkerFaceColor',gammaC{y},'MarkerEdgeColor',gammaC{y},'Color',gammaC{y})
plot([0 400],[0 400],'--')
[rr,pp] = corrcoef(lat_p100(:,1),lat_p100(:,2));
title('P100 peak time')
ylim([0 400])
xlim([0 400])
xlabel(sprintf('r = %2.2f, p = %0.2g',[rr(1,2) pp(1,2)]))


%% Bland-Altman plot

% Peak analysis
peak_diff_sess = diff(peak_diff,[],2);
peak_mean_sess = mean(peak_diff,2);

figure(302)
subplot(1,3,1)
hold on
plot(peak_mean_sess,peak_diff_sess,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5])
plot([min(peak_mean_sess) max(peak_mean_sess)],[mean(peak_diff_sess) mean(peak_diff_sess)],'--k')
plot([min(peak_mean_sess) max(peak_mean_sess)],[std(peak_diff_sess)*1.96 std(peak_diff_sess)*1.96],'--k')
plot([min(peak_mean_sess) max(peak_mean_sess)],[std(peak_diff_sess)*-1.96 std(peak_diff_sess)*-1.96],'--k')
title('Bland Altman N75-P100 peak difference')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(peak_mean_sess), max(peak_mean_sess)]; ax.YLim = [min(peak_diff_sess), max(peak_diff_sess)];


diff_lat_n75_sess = diff(lat_n75,[],2);
mean_lat_n75_sess = mean(lat_n75,2);

figure(302)
subplot(1,3,2)
hold on
plot(mean_lat_n75_sess,diff_lat_n75_sess,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5])
plot([min(mean_lat_n75_sess) max(mean_lat_n75_sess)],[mean(diff_lat_n75_sess) mean(diff_lat_n75_sess)],'--k')
plot([min(mean_lat_n75_sess) max(mean_lat_n75_sess)],[std(diff_lat_n75_sess)*1.96 std(diff_lat_n75_sess)*1.96],'--k')
plot([min(mean_lat_n75_sess) max(mean_lat_n75_sess)],[std(diff_lat_n75_sess)*-1.96 std(diff_lat_n75_sess)*-1.96],'--k')
title('Bland Altman P100 peak time')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_lat_n75_sess), max(mean_lat_n75_sess)]; ax.YLim = [min(diff_lat_n75_sess), max(diff_lat_n75_sess)];

diff_lat_p100_sess = diff(lat_p100,[],2);
mean_lat_p100_sess = mean(lat_p100,2);

figure(302)
subplot(1,3,3)
hold on
plot(mean_lat_p100_sess,diff_lat_p100_sess,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5])
plot([min(mean_lat_p100_sess) max(mean_lat_p100_sess)],[mean(diff_lat_p100_sess) mean(diff_lat_p100_sess)],'--k')
plot([min(mean_lat_p100_sess) max(mean_lat_p100_sess)],[std(diff_lat_p100_sess)*1.96 std(diff_lat_p100_sess)*1.96],'--k')
plot([min(mean_lat_p100_sess) max(mean_lat_p100_sess)],[std(diff_lat_p100_sess)*-1.96 std(diff_lat_p100_sess)*-1.96],'--k')
title('Bland Altman P100 peak time')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_lat_p100_sess), max(mean_lat_p100_sess)]; ax.YLim = [min(diff_lat_p100_sess), max(diff_lat_p100_sess)];



% Informed basis set
diff_n75bw_sess = diff(squeeze(bandwidth(:,:,1)),[],2);
mean_n75bw_sess = mean(squeeze(bandwidth(:,:,1)),2);

figure(302)
subplot(3,3,1)
hold on
plot(mean_n75bw_sess,diff_n75bw_sess,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5])
plot([min(mean_n75bw_sess) max(mean_n75bw_sess)],[mean(diff_n75bw_sess) mean(diff_n75bw_sess)],'--k')
plot([min(mean_n75bw_sess) max(mean_n75bw_sess)],[std(diff_n75bw_sess)*1.96 std(diff_n75bw_sess)*1.96],'--k')
plot([min(mean_n75bw_sess) max(mean_n75bw_sess)],[std(diff_n75bw_sess)*-1.96 std(diff_n75bw_sess)*-1.96],'--k')
title('Bland Altman N75-P100 peak difference')
xlabel('Average of sessions')
ylabel('Difference between sessions')
ax=gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [min(mean_n75bw_sess), max(mean_n75bw_sess)]; ax.YLim = [min(diff_n75bw_sess), max(diff_n75bw_sess)];
