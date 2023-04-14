% Adapted from Glover (1999) for modeling bold fMRI data, applied to VEP
% Developed by Carlyn Patterson Gentile, last edited November 11, 2022
% Inputs
% - t is time in ms
% - p contains the model parameters as follow:
% - p(1) is n1 that shifts the peak of the first gamma function
% - p(2) is t1/n1, which determines first the peak location (in ms)
% - p(3) is a1 that scales the first (negative) gamma function
% - p(4) is n2 that shifts the peak of the second gamma function
% - p(5) is t2/n2, which determines the second peak location (in ms)
% - p(6) is a2 that scales the second (positive) gamma function
% - p(7) is n3 that shifts the peak of the third gamma function
% - p(8) is t3/n3, which determines the third peak location (in ms)
% - p(9) is a3 that scales the third (negative) gamma function
% - p(10) is n4 that shifts the peak of the fourth gamma function
% - p(11) is t4/n4, which determines the fourth peak location (in ms)
% - p(12) is a4 that scales the fourth (positive) gamma function

% Outputs:
% -  vep_fit is the fitted data
% -----------------------------------------------------------------------
% Example of how to fit to EEG/VEP data: myFx compares model data to real
% data, p0 is model parameter guesses, yFit is to get the fitted data
% myFx = @(p) sqrt(sum((voltage - gammaVEP_model3(time,p)).^2));
% p0 = [35 75 min(voltage) 25 100 max(voltage) 27 135 min(voltage) 30 220 max(voltage)];
% Mdl_param = fmincon(myFx,p0,[],[],[],[],lb,ub);
% yFit = gammaVEP_model3(time,Mdl_param);

function [vep_fit,gamma1,gamma2,gamma3,gamma4,gamma5] = gammaVEP_model5(t,p)

% n1/t1 ~= 75, to capture the n75 peak
n1 = p(1);
t1 = p(2)/n1;
a1 = p(3);
c1 = 1/max((t.^n1).*exp(-t./t1));

% n2/t2 ~= 100, to capture the p100 peak
n2 = p(4);
t2 = p(5)/n2;
a2 = p(6);
c2 = 1/max((t.^n2).*exp(-t./t2));

% n3/t3 ~= 135, to capture the n135 peak
n3 = p(7);
t3 = p(8)/n3;
a3 = p(9);
c3 = 1/max((t.^n3).*exp(-t./t3));

% n4/t4 ~= 220, pos peak
n4 = p(10);
t4 = p(11)/n4;
a4 = p(12);
c4 = 1/max((t.^n4).*exp(-t./t4));

% n5/t5 ~= 300, neg peak
n5 = p(13);
t5 = p(14)/n5;
a5 = p(15);
c5 = 1/max((t.^n5).*exp(-t./t5));



gamma1 = zeros(1,length(t));
gamma2 = zeros(1,length(t));
gamma3 = zeros(1,length(t));
gamma4 = zeros(1,length(t));
gamma5 = zeros(1,length(t));

    for i = 1:length(t)
        gamma1(:,i) = (c1*(t(i)^n1)*(exp(-t(i)/t1)));
        gamma2(:,i) = (c2*(t(i)^n2)*(exp(-t(i)/t2)));
        gamma3(:,i) = (c3*(t(i)^n3)*(exp(-t(i)/t3)));
        gamma4(:,i) = (c4*(t(i)^n4)*(exp(-t(i)/t4)));
        gamma5(:,i) = (c5*(t(i)^n5)*(exp(-t(i)/t5)));
    end

    gamma1 = a1.*(gamma1./max(gamma1));
    gamma2 = a2.*(gamma2./max(gamma2));
    gamma3 = a3.*(gamma3./max(gamma3));
    gamma4 = a4.*(gamma4./max(gamma4));
    gamma5 = a5.*(gamma5./max(gamma5));

    vep_fit = gamma1 + gamma2 + gamma3 + gamma4 + gamma5;
 
end