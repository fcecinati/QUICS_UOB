s = cell(12,1);
m=11;
y=2007;
for i=1:12
    s(i)=cellstr([num2str(y),'-',num2str(m,'%02d'),'-01 00:00:00']);
    m=m+1;
    if m==13
        m=1;
        y=y+1;
    end
end
for i=1:12
    NewGeneratorRT( char(s(i)))
end

EnsemblesHBV

