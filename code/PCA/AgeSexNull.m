function [adjustedPCAscores]=AgeSexNull(nDims,PCAscore,sex,age)
% Function adjusts PCA scores to a sex nulled subject at age 15
% INPUT
% nDims - number of PCs to analyze
% PCAscore - the N-by-p matrix of PCA scores (N is subjects, p is PCs)
% sex - column vector identifying the sex of each subject (female = 1, male = 2)
% age - column vector identifying the age of each subject (in years)
% If you only want to adjust by sex or age, leave the other vector empty []
% Define synthetic age range

    for i=1:nDims
        if isempty(age)==0
            %Adjust each score by age
            temp=polyfit(age,PCAscore(:,i),1);
            adjusted_age=15;
            adjustAge_allSub=polyval(temp,adjusted_age)-polyval(temp,age);
        else
            adjustAge_allSub=0;
        end
        
        if isempty(sex)==0
            %Adjust each score by sex
            temp2=polyfit(sex,PCAscore(:,i),1);
            adjusted_sex=1.5; % this is half way between female (1) and male (2)
            adjustSex_allSub=polyval(temp2,adjusted_sex)-polyval(temp2,sex);
        else
            adjustSex_allSub=0;
        end

        %adjust score by age and sex
        adjustedPCAscores(:,i)=PCAscore(:,i)+adjustAge_allSub+adjustSex_allSub;
    end

end