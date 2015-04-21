mu = 0;
sigma = 5.6855;
h = (1:128)';
var_temp = zeros(10,10,10);
var_matrix = zeros(10,10);
i=0; j=0;

for i=1:10
    a=i*10;
    sigma1 = (sigma)/(1.09641257*(a^(-0.94003379)));
    for j=1:10
        b=1/(j*10);
        for k=1:10               
            %Generate normal random fields
            d=normrnd(mu,sigma,256,256);
            % Apply a lowpass filter
            base = 1:a;
            fun = ((sqrt(2/pi)*sin(pi/8)*gamma(3/4))./base.^(3/4)).^b;
            f1 = [fliplr(fun),1,fun];
            f2 = ftrans2(f1);
            dnew = filter2(f2,d)./a;

            var_temp(i,j,k) = var(dnew(:)); 
        end
        var_matrix = mean(var_temp,3);
    end
    
end
sigma_new = sqrt(var_matrix);
sig_matrix = sigma_new./sigma;
var_matrix = var_matrix./(sigma^2);
count=1:10;
aplot = count*10;
bplot = count*0.002;
linecolours = colormap('jet');

fig1 = figure;
figure(fig1);
hold on
for i=1:10
    plot(bplot,sig_matrix(i,:),'color',linecolours(i*6,:))
end
%ylim([0 2])
title('STD as function of b')
saveas(fig1,'\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test1\Sigma_fb_avg_corrected.jpg')

fig2 = figure;
figure(fig2);
hold on
for i=1:10
    plot(aplot, sig_matrix(:,i),'color',linecolours(i*6,:))
end
title('STD as function of a')
%ylim([0 2])
saveas(fig2,'\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test1\Sigma_fa_avg_corrected.jpg')




