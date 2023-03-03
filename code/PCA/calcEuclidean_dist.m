function [Euclidean_dist]=calcEuclidean_dist(PC_no,X1,X2)


for i=1:size(X1,1)
    for j=1:size(X2,1)
        for k=1:PC_no
            temp(i,j,k)=(X1(i,k)-X2(j,k))^2;
        end
    end
end

Euclidean_dist=sqrt(squeeze(sum(temp,3)));

end