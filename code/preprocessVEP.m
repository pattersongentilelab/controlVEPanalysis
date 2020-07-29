% Preprocess VEP for controls and cases

load VEP_and_subject_data vep_files vep

%% Organize and filter VEP data, and remove high amplitude trials

Fs=1024; % Sampling rate
reversal=2; % number of reversals

discarded_trials=0;
total_trials=0;
 
for x=1:size(vep,1)
    all_VEP=cell2mat(vep(x,3));
    
    duration=size(all_VEP,1)./Fs; % total duration of recording (in seconds)

    % Bandstop filter for 60Hz noise in VEP signal
    d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',58,...
        'HalfPowerFrequency2',62,'SampleRate',Fs);
    all_VEP=filter(d,all_VEP);

    d=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',119,...
        'HalfPowerFrequency2',121,'SampleRate',Fs);
    all_VEP=filter(d,all_VEP);
    
    % Break up data into individual reversals
    x_all_VEP=1/Fs:1/Fs:duration;
    x_data=1/Fs:1/Fs:1/reversal;

    y_data=reshape(all_VEP,[Fs/reversal,duration*reversal]);

    % remove trials where the absolute max voltage >1 or < 0.05
    temp=find(max(abs(y_data),[],1)<1 & max(abs(y_data),[],1)>0.05);
    bad_trials=size(y_data,2)-length(temp);
    discarded_trials=discarded_trials+bad_trials;
    total_trials=total_trials+size(y_data,2);
    disp(['number of discarded trials =' num2str(bad_trials)]);
    y_data=y_data(:,temp)';
    
    vep{x,3}=y_data;
    vep{x,4}=x_data;
    
    clear x_data y_data
end

disp(['proportion of discarded trials =' num2str(discarded_trials./total_trials)]);


%% Organize control subject VEP
control_loc=find(vep_files.subjecttype=='Control');

control_vep_temp=vep(control_loc,:);
control_vep_file_temp=vep_files(control_loc,:);

counter1=1;
counter2=1;
no_dup=unique(cell2mat(control_vep_temp(:,1)));

for x=1:size(no_dup,1)
    temp_ID=cell2mat(control_vep_temp(counter2,1));
    temp_loc=find(cell2mat(control_vep_temp(:,1))==temp_ID);
    temp_x=cell2mat(control_vep_temp(counter2,4));
    temp_vep=[];
    for y=1:length(temp_loc)
        temp_vep2=cell2mat(vep(temp_loc(y),3));
        temp_vep=cat(1,temp_vep,temp_vep2);
        clear temp_vep2
    end
    control_vep{counter1,1}=temp_ID;
    control_vep{counter1,2}=temp_vep;
    control_vep{counter1,3}=temp_x;
    control_vep_files(counter1,:)=control_vep_file_temp(temp_loc(1),:);
    counter1=counter1+1;
    counter2=counter2+length(temp_loc);
    clear temp*
end

% Normalize control data
for x=1:length(control_vep)
    y_data=cell2mat(control_vep(x,2));
    
    % Normalize by max voltage
    max_y=max(nanmedian(y_data,1));
    y_data=y_data./max_y;
    
    % Normalize pre-response to 0
    for y=1:size(y_data,1)
        temp=median(median(y_data(y,1:51)));
        y_data(y,:)=y_data(y,:)-(temp*ones(size(y_data(y,:))));
    end
    
    control_vep{x,2}=y_data;
end


clear *temp

%% Organize case subject VEP
case_loc=find(vep_files.subjecttype=='Case');

case_vep_temp=vep(case_loc,:);
case_vep_file_temp=vep_files(case_loc,:);

counter1=1;
counter2=1;
no_dup=unique(cell2mat(case_vep_temp(:,1:2)),'rows');

for x=1:size(no_dup,1)
    temp_ID=cell2mat(vep(counter2,1));
    temp_dayspostinj=cell2mat(case_vep_temp(counter2,2));
    temp_loc=find(cell2mat(case_vep_temp(:,1))==temp_ID & cell2mat(case_vep_temp(:,2))==temp_dayspostinj);
    temp_x=cell2mat(case_vep_temp(counter2,4));
    temp_vep=[];
    for y=1:length(temp_loc)
        temp_vep2=cell2mat(vep(temp_loc(y),3));
        temp_vep=cat(1,temp_vep,temp_vep2);
        clear temp_vep2
    end
    case_vep{counter1,1}=temp_ID;
    case_vep{counter1,2}=temp_dayspostinj;
    case_vep{counter1,3}=temp_vep;
    case_vep{counter1,4}=temp_x;
    case_vep_files(counter1,:)=case_vep_file_temp(temp_loc(1),:);
    counter1=counter1+1;
    counter2=counter2+length(temp_loc);
    clear temp*
end

% Normalize case data
for x=1:length(case_vep)
    y_data=cell2mat(case_vep(x,3));
    
    % Normalize by max voltage
    max_y=max(nanmedian(y_data,1));
    y_data=y_data./max_y;
    
    % Normalize pre-response to 0
    for y=1:size(y_data,1)
        temp=median(median(y_data(y,1:51)));
        y_data(y,:)=y_data(y,:)-(temp*ones(size(y_data(y,:))));
    end
    
    case_vep{x,3}=y_data;
end
    
clear *temp

save cleaned_VEP control_vep control_vep_files case_vep case_vep_files
