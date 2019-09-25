%coin flip

clear all close all;
t=0;
for u = [1 10 100 1000 10000 100000 1000000];
    t = t+1;
    j=0;
    for i = 1:u;
        j = j+1;
        samples = 100;
        a= round(rand(samples,1));
        b=round(rand(samples,1));
        c = a == b;
        d = sum(c)/samples;
        e(j) = d;
    end
    f(t)= std(e);
end

plot(f);


%%
n = 1+14.88*((.05/.01)^2)