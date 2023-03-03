% function to calculate the width of 1/2 height of a gamma function
% input
% - minT is the start time of the gamma function (in ms)
% - maxT is the end time of the gamma function (in ms)
% - mdl_fit is a 1x2 matrix containing n1(1), n1*t1(2) for the
% gamma function

function [bandwidth] = gamma_bandwidth(minT,maxT,mdl_fit)

t = minT:0.01:maxT;
n1 = mdl_fit(1);
t1 = mdl_fit(2)/mdl_fit(1);

gamma_fx = ((t.^n1).*(exp(-t./t1)));

if isinf(max(gamma_fx)) == 1
    gamma_fx = gamma_fx(~isinf(gamma_fx));
    t = t(~isinf(gamma_fx));
end

if max(gamma_fx)<=0
    gamma_fx = gamma_fx.*-1;
end

gamma_fx = gamma_fx./max(gamma_fx);


Tgamma50 = t(gamma_fx>=0.5);

bandwidth = Tgamma50(end) - Tgamma50(1);

end