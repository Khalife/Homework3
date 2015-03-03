function [X,accpt] = nsrwHM(x0,N,c,d,Gamma)
% function [X,alpha] = nsrwHM(x0,N)
%   X the HM sequence simulating Pi
%   naive symetric random walk
%   x0 the initial d-dim. value
%   N the number of iterations
    
    u=0; % this variable will store the uniform distribution to simuate a bernoulli dist.
    accpt=zeros(N,1);%the acceptation rate
    alpha=zeros(N,1); 
    Y=zeros(1,d); % the proposal
    
    X=zeros(N,d);
    X(1,:)=x0;
    
    S=diag(Gamma)';
    Z=zeros(1,d);
    
    for n=2:N
        Y=X(n-1,:)+c*normrnd(0,1,1,d);
        aux=normpdf(Y,Z,S)/normpdf(X(n-1,:),Z,S);
        alpha(n)=min([1,aux]);
        
        u=rand;
        if u<alpha(n)
            X(n,:)=Y;
            accpt(n)=1;
        else
            X(n,:)=X(n-1,:);
            accpt(n)=0;
        end
    end
    
end