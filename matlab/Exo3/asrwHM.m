function [X,accpt,lambda,subopt] = asrwHM(x0,N,n0,c,d,Gamma)
% function [X,alpha] = nsrwHM(x0,N)
%   X the HM sequence simulating Pi
%   adaptive symetric random walk
%   x0 the initial d-dim. value
%   N the number of iterations
    
    u=0; % this variable will store the uniform distribution to simuate a bernoulli dist.
    accpt=zeros(N,1);%the acceptation rate
    alpha=zeros(N,1); 
    Y=zeros(1,d); % the proposal
    
    lambda=zeros(N-n0,d);
    subopt=zeros(N-n0,1);
    
    X=zeros(N,d);
    X(1,:)=x0;
    
    S=diag(Gamma)';
    Z=zeros(1,d);
    
    C=(2.38)^2;
    
    for n=2:n0
        Y=X(n-1,:)+c*normrnd(0,1,1,d);
        aux=normpdf(Y,Z,S)/normpdf(X(n-1,:),Z,S);
        alpha(n)=min([1,aux]);
        
        u=rand;
        if u<alpha(n)
            X(n,:)=Y;
            accpt(n)=0;
        else
            X(n,:)=X(n-1,:);
            accpt(n)=1;
        end
    end
    
    m=1/n0*sum(X(1:n0));
    EGamma=(X(n0,:)-m)'*(X(n0,:)-m)./n0;
    
    for n=n0+1:N
        Y=X(n-1,:)+mvnrnd(Z,C/d*EGamma);
        aux=mvnpdf(Y,Z,S)/mvnpdf(X(n-1,:),Z,S);
        alpha(n)=min([1,aux]);
        
        u=rand;
        if u<alpha(n)
            X(n,:)=Y;
            accpt(n)=1;
        else
            X(n,:)=X(n-1,:);
            accpt(n)=0;
        end
        
        m=m+1/(n-n0+1).*(X(n,:)-m);
        EGamma=((X(n0,:)-m)'*(X(n0,:)-m)-EGamma)./(n-n0+1)+EGamma;
        
        [Q,D]=eig(EGamma);
        D=diag(sqrt(abs(diag(D))));
        A=Q*D*Q';
        
        M=A/sqrtm(Gamma);
        lambda(n-n0,:)=eig(M);
        denom=sum(lambda(n-n0,:).^-2);
        num=sum(lambda(n-n0,:).^-1)^2;
        subopt(n-n0)=d*denom/num;
        
    end
    
end