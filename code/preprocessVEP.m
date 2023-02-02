% Preprocess VEP for controls and cases

data_path = getpref('controlVEPanalysis','MindsMatter_DataPath');

load([data_path '/VEP_and_subject_data.mat'],'vep_files','raw_vep')

%% Organize and filter VEP data, and remove high and low amplitude trials

Fs=1024; % Sampling rate
reversal=2; % number of reversals

discarded_trials=0;
total_trials=0;
 
for x=1:size(raw_vep,1)
    all_VEP=cell2mat(raw_vep(x,3));
    
    duration=size(all_VEP,1)./Fs; % total duration of recording (in seconds)

    % Bandstop filter for 60Hz noise in VEP signal
    d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',59,...
        'HalfPowerFrequency2',61,'SampleRate',Fs);
    all_VEP=filter(d,all_VEP);

    d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',119,...
        'HalfPowerFrequency2',121,'SampleRate',Fs);
    all_VEP=filter(d,all_VEP);
    
    % Break up data into individual reversals
    x_all_VEP=1/Fs:1/Fs:duration;
    x_data=1/Fs:1/Fs:1/reversal;

    y_data=reshape(all_VEP,[Fs/reversal,duration*reversal]);

    % remove trials where the absolute max voltage >1
    temp=find(max(abs(y_data),[],1)<1);
    temp2=find(max(abs(y_data),[],1)>=1);
%     figure(10)
%     plot(y_data(:,temp2))
%     title('bad trial')
%     pause
    bad_trials=size(y_data,2)-length(temp);
    discarded_trials=discarded_trials+bad_trials;
    total_trials=total_trials+size(y_data,2);
    disp(['number of discarded trials =' num2str(bad_trials)]);
    y_data=y_data(:,temp)';
    
    VEP{x,1}=cell2mat(raw_vep(x,1));
    VEP{x,2}=cell2mat(raw_vep(x,2));
    VEP{x,3}=x_data;
    VEP{x,4}=y_data;
    
    clear x_data y_data
end

disp(['proportion of discarded trials =' num2str(discarded_trials./total_trials)]);


%% Organize VEP data into individual sessions

counter1=1;
counter2=1;
no_dup=unique(cell2mat(VEP(:,1:2)),'rows');

for x=1:size(no_dup,1)
    temp_ID=no_dup(x,1);
    temp_age=no_dup(x,2);
    temp_loc=find(cell2mat(VEP(:,1))==temp_ID & cell2mat(VEP(:,2))==temp_age);
    temp_loc2=find(table2array(vep_files(:,1))==temp_ID & table2array(vep_files(:,212))==temp_age);
    cleaned_vep_files_loc(x,:)=temp_loc2(1);
    temp_x=cell2mat(VEP(counter2,3));
    temp_vep=[];
    for y=1:length(temp_loc)
        temp_vep2=cell2mat(VEP(temp_loc(y),4));
        temp_vep=cat(1,temp_vep,temp_vep2);
        clear temp_vep2
    end
    cleaned_vep{counter1,1}=temp_ID;
    cleaned_vep{counter1,2}=temp_age;
    cleaned_vep{counter1,3}=temp_x;
    cleaned_vep{counter1,4}=temp_vep;
    counter1=counter1+1;
    counter2=counter2+length(temp_loc);
    clear temp_loc temp_ID temp_age temp_vep
end

cleaned_vep_files=vep_files(cleaned_vep_files_loc,:);

% Normalize data
for x=1:length(cleaned_vep)
    y_data=cell2mat(cleaned_vep(x,4));
    
%     % Normalize by max voltage
%     max_y=max(nanmedian(y_data,1));
%     y_data=y_data./max_y;
    
    % Normalize pre-response to 0
    for y=1:size(y_data,1)
        temp=mean(mean(y_data(y,1:51)));
        y_data(y,:)=y_data(y,:)-(temp*ones(size(y_data(y,:))));
    end
    
    cleaned_vep{x,4}=y_data;
end


clear *temp

% save cleaned_VEP cleaned_vep cleaned_vep_files

clear
