% Adapted from Glover (1999) for modeling bold fMRI data, applied to VEP
% Developed by Carlyn Patterson Gentile, last edited July 13, 2023
% Inputs
% - t is time in ms
% - p contains the model parameters in groups of 3 as follows:
% - p(1,4,6,8) is n1 that shifts the peak of the gamma function and
% impacts bandwidth
% - p(2) is t1/n1 for gamma 1 only, which determines the peak location (in ms)
% - p(3,5,7,9) is a1 that scales the gamma function
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

function [vep_fit,gamma] = gammaVEP_modelReduced(t,p)

nGamma = 4;
gamma = zeros(nGamma,length(t));

for x = 1:nGamma
    % get the parameters for gamma x
    switch x
        case 1
            a1 = p(3);
            n1 = p(1);
            t1 = p(2)/n1;
        case 2
            a1 = p(5);
            n1 = p(4);
            t1 = -0.4*(nPrec) + 126;
        case 3
            a1 = p(7);
            n1 = p(6);
            t1 = -0.8*(nPrec) + 186;
        case 4
            a1 = p(9);
            n1 = p(8);
            t1 = -1.1*(nPrec) + 260;
    end
    
    c1 = 1/max((t.^n1).*exp(-t./t1));


    for i = 1:length(t)
        gamma(x,i) = (c1*(t(i)^n1)*(exp(-t(i)/t1)));
    end

    gamma(x,:) = a1.*(gamma(x,:)./max(gamma(x,:)));
    nPrec = n1;
end

vep_fit = sum(gamma,1);
end