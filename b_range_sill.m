% relationship between b and the variog parameters
clear all

tic
btest = 30:10:150; bs = size(btest,2);
atest = 20:10:200; as = size(atest,2);
stdtest = 1:10;  stds = size(stdtest,2);
di=1:128;  ds = size(di,2);

variog = zeros(ds,1);
s=0; c=0; j=0;
sill = zeros(as,bs,stds,10);
range = zeros(as,bs,stds,10);
nugget = zeros(as,bs,stds,10);

for stdev = stdtest
    for i = 1:10
        d = normrnd(0,stdev,256,256);
        s=s+1; c=0;
        for ai = atest
            c=c+1; j=0;
            for bi = btest
                j=j+1;            
                base = 1:ai;
                fun = exp(-(base./bi));
                f1 = [fliplr(fun),(1),fun];
                H0 = sum(f1);
                f1=f1/H0;
                f2 = ftrans2(f1);
                dnew = filter2(f2,d);
                newstd = std(dnew(:));
                dnew = dnew*stdev/newstd;

                [hvar, vvar] = Variogram(dnew);
                variog = (hvar + vvar)/2;
                myfunction2 = fittype('c.*(1-exp(-(h./a)))+b','independent',{'h'},'coefficients',{'c','b','a'});
                myoption2 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[stdev^2 10 50]);
                try
                    [myfit2, gof] = fit(di',variog,myfunction2,myoption2);
                    sill(c,j,s,i) = myfit2.c;
                    range(c,j,s,i) = myfit2.a;
                    nugget(c,j,s,i) = myfit2.b;
                catch
                    sill(c,j,s,i) = NaN;
                    range(c,j,s,i) = NaN;
                    nugget(c,j,s,i) = NaN;
                end        
            end
        end
    end  
    sill = nanmean(sill,4);
    range = nanmean(sill,4);
    nugget = nanmean(nugget,4);
    s=0;
end

save('test','atest','btest','stdtest','range','sill','nugget')

f1=figure;
figure(f1)
set(f1,'Position',[200,200,600,300])
for i=1:10
    imagesc(atest,btest,range(:,:,i))
    xlabel('atest')
    ylabel('btest')
    colorbar
    name = ['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\CalibrationTest\range-std',num2str(i)];
    saveas(f1,name,'png')
end

f2=figure;
figure(f2)
set(f2,'Position',[200,200,600,300])
for i=1:10
    imagesc(atest,btest,sill(:,:,i))
    xlabel('atest')
    ylabel('btest')
    colorbar
    name = ['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\CalibrationTest\sill-std',num2str(i)];
    saveas(f2,name,'png')
end

f3=figure;
figure(f3)
set(f3,'Position',[200,200,600,300])
for i=1:10
    imagesc(atest,btest,nugget(:,:,i))
    xlabel('atest')
    ylabel('btest')
    colorbar
    name = ['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\CalibrationTest\nugget-std',num2str(i)];
    saveas(f3,name,'png')
end

