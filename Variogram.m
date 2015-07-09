function [variog] = Variogram(M) 

sz = size(M);

%Generation of random sample
x = ceil(rand(1000,1).*sz(1));
y = ceil(rand(1000,1).*sz(2));
sample = zeros(1000,1);

for i=1:1000
    sample(i) = M(x(i),y(i));
end

%Calculating distances
distance = zeros(1000);
for x1=1:1000
    for x2=1:1000
        distance(x1,x2) = (((x(x1)-x(x2))^2)+((y(x1)-y(x2))^2))^0.5;
    end
end
dis = ceil(max(distance(:)));

for i=2:dis   
    [c,l]=find(distance>=(i-1)&distance<i);
    mat = ((sample(c)-sample(l)).^2);
    variog(i-1)= nanmean(mat(:));
end
