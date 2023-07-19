% Adapted from Glover (1999) for modeling bold fMRI data, applied to VEP
% Developed by Carlyn Patterson Gentile, last edited July 13, 2023
% Inputs
% - t is time in ms
% - p contains the model parameters in groups of 3 as follows:
% - p(1,4,7...) is n1 that shifts the peak of the gamma function
% - p(2,5,8...) is t1/n1, which determines the peak location (in ms)
% - p(3,6,9...) is a1 that scales the gamma function
% - nGamma is the number of gamma functions being included in the model
% Outputs:
% -  vep_fit is the fitted data
% -  Gamma holds the individual gamma functions
% -----------------------------------------------------------------------
% Example of how to fit to EEG/VEP data: myFx compares model data to real
% data, p0 is model parameter guesses, yFit is to get the fitted data
% myFx = @(p) sqrt(sum((voltage - gammaVEP_model3(time,p)).^2));
% p0 = [35 75 min(voltage) 25 100 max(voltage) 27 135 min(voltage) 30 220 max(voltage)];
% Mdl_param = fmincon(myFx,p0,[],[],[],[],lb,ub);
% yFit = gammaVEP_model3(time,Mdl_param);

function [vep_fit,gamma] = gammaVEP_model(t,p,nGamma)

gamma = zeros(nGamma,length(t));
y = 1;

for x = 1:nGamma
    % get the parameters for gamma x
    n1 = p(y);
    t1 = p(y+1)/n1;
    a1 = p(y+2);
    c1 = 1/max((t.^n1).*exp(-t./t1));


    for i = 1:length(t)
        gamma(x,i) = (c1*(t(i)^n1)*(exp(-t(i)/t1)));
    end

    gamma(x,:) = a1.*(gamma(x,:)./max(gamma(x,:)));
    y = y+3;
end

vep_fit = sum(gamma,1);
end