function var_explained=calcVarExplained(ydata,coeff,PCAmodel)


for y=1:size(PCAmodel,2)
    temp_mdl=PCAmodel(:,y);
    temp_var=var(temp_data-temp_mdl');
    temp_var_diff=tempTvar-temp_var;
    tempTvar_mdl=tempTvar_mdl-temp_var_diff;
    var_unexplained_subject(x,y)=tempTvar_mdl./tempTvar;
end

end