% Exercise 1 : gamma Metropolis hastings

%try different values of gamma
n_samples=10000;
%gamma=1;
gamma=4;
%gamma=50;
%[error,X_sample]=gammaMH(1,n_samples);
mcmc=10;
%loop to see the influence of gamma
Gamma=[1;4;50];
error=zeros(3,1);

for i=1:3
    for n=1:mcmc
        [e,~]=gammaMH(Gamma(i),n_samples);
        error(i)=error(i)+e;
        n
    end
    i
end
error=error/mcmc;
plot(error)
title('Error plot');