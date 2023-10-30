
function [tbl] = lmfit_tbl(mdl)

predicted = mdl.Coefficients.Estimate(2:end);
se = mdl.Coefficients.SE(2:end);
p_val = mdl.Coefficients.pValue(2:end);
VarNames = mdl.Coefficients.Properties.RowNames(2:end);
tbl = table('Size',[length(predicted) 1],'VariableTypes',{'double'},'VariableNames',{'estimate'},'RowNames',VarNames);

tbl.estimate = predicted;
tbl.low95 = predicted - (1.96*se);
tbl.hi95 = predicted + (1.96*se);
tbl.p_val = p_val;
