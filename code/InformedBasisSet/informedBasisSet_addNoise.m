% load train and test subjects 2021 TVST PCA paper, 5 gamma model
analysis_path = getpref('controlVEPanalysis','controlVEP_AnalysisPath');
load([analysis_path '/controlTrainTest'])

xdata = xdata.*1000; % scale time to ms

Vep = [squeeze(nanmean(control_train_vep,2));squeeze(nanmean(control_test_vep,2))].*100; % x100 to correct scaling error on diopsys voltage output
subject_data = [control_train;control_test];


for I = 1:4
    
        % introduce noise to determine the stability of the gamma functions
    switch I
        case 1 % raw data
            vep = Vep;
            plotColor = 'k';
            plot_shift = -0.3;
        case 2 %random noise
                r_noise = rand(size(vep)).*10;
                vep = Vep + r_noise;
                plotColor = 'r';
                plot_shift = -0.1;
        case 3 % waveform shift
            vep = cat(2,Vep(:,20:end),Vep(:,1:19)); % shifts the waveform
            plotColor = 'b';
            plot_shift = 0.1;
        case 4 % remove bounds
            vep = Vep;
            plotColor = 'm';
            plot_shift = 0.3;  
    end
   

    % truncate the part of the VEP fit to the model to improve fits
    minF = 1;
    maxF = 400;
    vepF = vep(:,minF:maxF);
    xdataF = xdata(minF:maxF);

    % Determine fits on individual VEP data
    yFit = zeros(size(vepF));
    gamma1 = zeros(size(vepF));
    gamma2 = zeros(size(vepF));
    gamma3 = zeros(size(vepF));
    gamma4 = zeros(size(vepF));
    gamma5 = zeros(size(vepF));
    mdl = zeros(size(vepF,1),15);
    r2 = zeros(size(vepF,1),1);
    bandwidth1 = zeros(size(vepF,1),1);
    bandwidth2 = zeros(size(vepF,1),1);
    bandwidth3 = zeros(size(vepF,1),1);
    bandwidth4 = zeros(size(vepF,1),1);
    bandwidth5 = zeros(size(vepF,1),1);

    figure
    for i = 1:size(vepF,1)
        ydata = vepF(i,:);
        myFx = @(p) sqrt(sum((ydata - gammaVEP_model5(xdataF,p)).^2));
        % Model guess
        switch I
            case 4
                p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 0 0 0];
                lb = [0 0 -Inf 0 0 0 0 0 -Inf 0 0 0 0 0 -Inf]; 
                ub = [Inf Inf 0 Inf Inf Inf Inf Inf 0 Inf Inf Inf 0.01 Inf 0];
            otherwise
                p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata)];
                lb = [5 50 -50 5 70 0 5 90 -50 5 200 0 5 300 -50]; 
                ub = [100 110 0 100 150 50 100 200 0 100 350 50 100 500 0];
        end
        mdl(i,:) = fmincon(myFx,p0,[],[],[],[],lb,ub);
        [vep_fit] = gammaVEP_model5(0:0.1:500,mdl(i,:));
        bandwidth1(i) = gamma_bandwidth(0,500,mdl(i,1:2));
        bandwidth2(i) = gamma_bandwidth(0,500,mdl(i,4:5));
        bandwidth3(i) = gamma_bandwidth(0,500,mdl(i,7:8));
        bandwidth4(i) = gamma_bandwidth(0,500,mdl(i,10:11));
        bandwidth5(i) = gamma_bandwidth(0,500,mdl(i,13:14));
        [yFit(i,:),gamma1(i,:),gamma2(i,:),gamma3(i,:),gamma4(i,:),gamma5(i,:)] = gammaVEP_model5(xdataF,mdl(i,:));
        r = corrcoef(ydata,yFit(i,:));
        r2(i,:) = r(1,2)^2;

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
        plot(xdataF,vepF(i,:),'.k')
        hold on
        plot(xdataF,yFit(i,:),['-' plotColor],'LineWidth',2)
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
        xlabel(sprintf('r2 = %2.2f',r2(i,:)))
        ylabel(sprintf('%2.0f',i))
        title(sprintf(['A1 = %2.0f, P1 = %2.1f, W1 = %2.0f; A2 = %2.0f, P2 = %2.1f, W2 = %2.0f;' ...
            'A3 = %2.0f, P3 = %2.1f, W3 = %2.0f; A4 = %2.0f, P4 = %2.1f, W4 = %2.0f; A5 = %2.0f, P5 = %2.1f, W5 = %2.0f'],...
            [mdl(i,3) mdl(i,2) bandwidth1(i) mdl(i,6) mdl(i,5) bandwidth2(i) mdl(i,9) mdl(i,8) bandwidth3(i) ...
            mdl(i,12) mdl(i,11) bandwidth4(i) mdl(i,15) mdl(i,14) bandwidth5(i)]))

        subplot(4,4,j+4)
        plot(xdataF,gamma1(i,:),'r')
        hold on
        plot(xdataF,gamma2(i,:),'b')
        plot(xdataF,gamma3(i,:),'g')
        plot(xdataF,gamma4(i,:),'m')
        plot(xdataF,gamma5(i,:),'c')
        ax=gca; ax.TickDir = 'out'; ax.Box = 'off';
    end

    % Determine fit on mean VEP data across individuals
    meanVEPf = mean(vepF,1);
    ydata = meanVEPf;
    % Model guess
    switch I
        case 4
            p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 0 0 0];
        otherwise
            p0 = [35 75 min(ydata) 25 100 max(ydata) 27 135 min(ydata) 30 220 max(ydata) 30 300 min(ydata)];
    end
    myFx = @(p) sqrt(sum((ydata - gammaVEP_model5(xdataF,p)).^2));
    Mdl = fmincon(myFx,p0,[],[],[],[],lb,ub);
    [yFit_m,gamma1m,gamma2m,gamma3m,gamma4m,gamma5] = gammaVEP_model5(xdataF,Mdl);
    figure
    plot(xdataF,meanVEPf,'.k')
    hold on
    plot(xdataF,yFit_m,['-' plotColor],'LineWidth',2)
    ax=gca; ax.TickDir = 'out'; ax.Box = 'off';

    bandwidth = cat(2,bandwidth1,bandwidth2,bandwidth3,bandwidth4,bandwidth5);

    % plotting parameters
    param = {'width 1','peak time 1','amplitude 1','width 2','peak time 2','amplitude 2','width 3','peak time 3','amplitude 3','width 4','peak time 4','amplitude 4','width 5','peak time 5','amplitude 5'};
    figure(200)
    y = 1;
    for x = 1:15
        subplot(5,3,x)
        switch x
            case [1,4,7,10]
                bw = bandwidth(:,y);
                hold on
                errorbar(1+plot_shift,mean(tDn),std(tDn),plotColor,'LineWidth',2,'MarkerFaceColor',plotColor,'MarkerEdgeColor',plotColor)
                y = y+1;
            otherwise
                hold on
                errorbar(1+plot_shift,mean(mdl(:,x)),std(mdl(:,x)),'o','LineWidth',2,'MarkerFaceColor',plotColor,'MarkerEdgeColor',plotColor)
        end
        ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1; ax.XTickLabel = param(x);
    end
    
    % Plot R2 values
    figure(202)
    hold on
    errorbar(1+plot_shift,mean(r2),std(r2),plotColor,'LineWidth',2,'MarkerFaceColor',plotColor,'MarkerEdgeColor',plotColor)

    % Plot params by age
    figure(201)
    for x = 1:15
        subplot(5,3,x)
        switch x
            case [1,4,7,10]
                tDn = mdl(:,x)./mdl(:,x-1);
                plot(subject_data.age_vep,tDn,['.' plotColor])
                lsline
            otherwise
                plot(subject_data.age_vep,mdl(:,x),['.' plotColor])
                lsline
        end
        ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; title(param(x));
    end
end
