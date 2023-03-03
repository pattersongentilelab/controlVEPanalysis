function [dprime p_val]=calcDprime(A,B)


%% calculate bootstrapped d-prime
Ams=bootstrp(1000,@(x)[mean(x) std(x)],A);
Bms=bootstrp(1000,@(x)[mean(x) std(x)],B);

for x=1:1000
    Astd2=Ams(x,2)^2;
    Bstd2=Bms(x,2)^2;
    ABstd2=(Astd2+Bstd2)/2;
    pooledSTD=sqrt(ABstd2);
    ABm=Ams(x,1)-Bms(x,1);
    bootD(x,:)=ABm/pooledSTD;
    clear pooledSTD Astd2 Bstd2 ABstd2 ABm
end

bootD=sort(bootD);
dprime=bootD([25 500 975],:);

clear ABm pooledSTD ABstd2

%% calculate d-prime null prediction
AB=cat(1,A,B);

ABms1=bootstrp(length(AB),@(x)[mean(x) std(x)],AB);
ABms2=bootstrp(length(AB),@(x)[mean(x) std(x)],AB);

for x=1:length(AB)
    ABstd2_1=ABms1(x,2)^2;
    ABstd2_2=ABms2(x,2)^2;
    ABstd2=(ABstd2_1+ABstd2_2)/2;
    pooledSTD=sqrt(ABstd2);
    ABm=ABms1(x,1)-ABms2(x,1);
    bootD_null(x,:)=ABm/pooledSTD;
    clear pooledSTD Astd2 Bstd2 ABstd2 ABm
end

bootD_null=sort(bootD_null);
dprime_null=bootD_null(round(length(AB)/2),:);

%% determine p-value
if dprime(2,:)<0
    dprime_null=-1*abs(dprime_null);
    p_val=length(find(bootD>dprime_null))./1000;
else
    dprime_null=abs(dprime_null);
    p_val=length(find(bootD<dprime_null))./1000;
end




end