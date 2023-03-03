function [simVEP]=synthVEP(nDims,PCAscore,PCAcoef,PCAweight,varargin)
% Function creates simulated VEP signals based on PCA for a 10, 15, and 20
% year old sex-nulled subject
% INPUT
% nDims - number of PCs to analyze
% PCAscore - the N-by-p matrix of PCA scores (N is subjects, p is PCs)
% PCAcoef - the m-by-p matrix for PCA coefficents (m is VEP (mV), p is PCs) 
% PCAweight - a row vector containing the variabilty accounted for by each PC


q = inputParser;
q.addParameter('age',[],@isnumeric); %age of subjects (years), if present will adjust VEP by age
q.addParameter('sex',[],@isnumeric); %sex of subjects, female = 1, male = 2, if present will adjust VEP by sex
q.addParameter('conc',[],@isnumeric); %concussion history, no history = 0, history = 1, if present will adjust VEP by concussion hx
q.parse(varargin{:});


% 
ages=[10 15 20];
sexes=[1 2];
conc_hx=[0 1];
C=['r';'b';'k'];

if isempty(q.Results.age)==0

    for i=1:length(ages)
        for j=1:nDims

            %Determine age adjustment
            temp=polyfit(q.Results.age,PCAscore(:,j),1);
            adjusted_age=ages(i);
            adjustAge=polyval(temp,adjusted_age);

            temp3=adjustAge.*PCAcoef(:,j)';
            simPC(j,:)=temp3;

            simVEP(i,:)=sum(simPC,1);
            x_data=1/1024:1/1024:0.5;
            
            if j==nDims
                figure(37)
                hold on
                plot(x_data,simVEP(i,:),'-','Color',C(i));
                ax=gca; ax.TickDir='out'; ax.Box='off';ax.XLim=[0 0.5];
                xlabel('Time (s)')
                ylabel('VEP signal (mV)')
                legend('10 years','15 years','20 years')
            end
        end
    end
end



if isempty(q.Results.sex)==0

    for i=1:length(sexes)
        for j=1:nDims

            %Determine sex adjustment
            temp=polyfit(q.Results.sex,PCAscore(:,j),1);
            adjusted_sex=sexes(i);
            adjustSex=polyval(temp,adjusted_sex);

            temp3=adjustSex*PCAcoef(:,j)'.*PCAweight(:,j)';
            simPC(j,:)=temp3;
        end
        simVEP(i,:)=sum(simPC,1);
        x_data=1/1024:1/1024:0.5;
        plot(x_data,simVEP(i,:),'-','Color',C(i));
        ax=gca; ax.TickDir='out'; ax.Box='off';ax.XLim=[0 0.5];ax.YLim=[-0.1 0.05];
        xlabel('Time (s)')
        ylabel('VEP signal (mV)')
        legend('female','male')
    end
end


if isempty(q.Results.conc)==0

    for i=1:length(conc_hx)
        for j=1:nDims

            %Determine sex adjustment
            temp=polyfit(q.Results.conc,PCAscore(:,j),1);
            adjusted_conc=conc_hx(i);
            adjustConc=polyval(temp,adjusted_conc);

            temp3=adjustConc*PCAcoef(:,j)';
            simPC(j,:)=temp3.*PCAweight(:,j)';

         simPC(j,:)=temp3;
        end
        simVEP(i,:)=sum(simPC,1);
        x_data=1/1024:1/1024:0.5;
        plot(x_data,simVEP(i,:),'-','Color',C(i));
        ax=gca; ax.TickDir='out'; ax.Box='off';ax.XLim=[0 0.5];ax.YLim=[-0.1 0.05];
        xlabel('Time (s)')
        ylabel('VEP signal (mV)')
        legend('no concussion hx','concussion hx')
    end
end



end