function [error,X_sample]=gammaMH(gamma,n_samples)
Pi_=[3/6 2/6 1/24 1/24 1/24 1/24];
Pi =@(x)(x > 0)*(x < 7)*Pi_(1,max(mod(x,7),1)); %extended distribution

X_set=zeros(1,gamma-1);
alpha=@(x,y)min(1,Pi(y)/Pi(x));

N=10*n_samples;
X=zeros(N,1);
Y=zeros(N,1);
x0=3;
X(1)=x0;

for n=2:N
    
%set of candidates for Y    
for g=0:(gamma-1)
    X_set(1,g+1)=X(n-1)-gamma+g;
end

for g=1:gamma
    X_set(1,g+gamma)=X(n-1)+g;
end

Y(n)=X_set(1,randi(2*gamma,1,1)); % y sampled with uniform distribution over the set
 
acceptation = rand(1,1);
if acceptation < alpha(X(n-1),Y(n)) 
X(n)=Y(n);
else     
X(n)=X(n-1);    
end

end

X_sample=X((N-n_samples):N,1); %to take in account burning period
values=[1 2 3 4 5 6];
theoretical_mean=sum(values.*Pi_);
theoretical_var=sum((values-theoretical_mean).^2.*Pi_);
error = abs(theoretical_var -var(X_sample));
%error2 = abs(theoretical_mean-mean(X_sample))
end
