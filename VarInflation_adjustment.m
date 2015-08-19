function [newens] = VarInflation_adjustment(raddata,ens)

sz = size(ens);
side = sz(1);
n_ens = sz(3);
newens = zeros(side, side, n_ens);

% Calculate discretization based on distance
x = (-(side/2)+0.5):((side/2)-0.5); y=x;
[X,Y] = meshgrid(x,y);
dist = zeros(side,side);
dist = (X.^2 + Y.^2).^0.5;
steps = 5:5:floor((side/2)*(2^0.5));
discret = zeros(side,side);
j=0;
for i=steps
    j=j+1;
    discret(dist>(i-5)&dist<i) = j;
end

% FFT radar data and ensembles
R = fft2(raddata); R = fftshift(R);
E = zeros(side,side,n_ens); 
for en = 1:n_ens
    ensemb = ens(:,:,en);
    E(:,:,en) = fft2(ensemb); E(:,:,en)=fftshift(E(:,:,en));
end

% Consider the mean of the ensembles FFT to calculate the ratio factors
Emean = mean(E,3);
correction = discret;
for i=1:j
    nom = abs(Emean(discret==i)); den = abs(R(discret==i));
    ratio(i) = mean(nom)/mean(den);
    correction(correction == i) = ratio(i);
end

% Apply the ratio factors to the ensembles FFT
for en = 1:n_ens
    Enew = E(:,:,en)./correction; 
    m = nanmin(Enew(:)); Enew(isnan(Enew)==1)=m; Enew(isinf(Enew)==1)=m;
    Enew = fftshift(Enew);
    newens(:,:,en) = ifft2(Enew,'symmetric');
end
newens(newens<0)=0;

        

