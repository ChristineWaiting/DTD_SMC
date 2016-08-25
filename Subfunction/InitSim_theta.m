function [theta]=InitSim_theta(Nsim,mu,sigma,dist,bound)
%
%---------------------------------------------------------
% Input  :   Nsim  : # of simulation
%            mu    : Nfirm*1
%            sigma : Nfirm*1
%            dist  : Nfirm*1
%            bound : Nfirm*[Low_bound up_bound] 
%--------------------------------------------------------- 
% Output : theta : Nsim* Nfirm
%---------------------------------------------------------
if nargin==4
    bound=[];
end
numP=length(mu);
theta = zeros(Nsim,numP);
RV=randn(Nsim,numP);
%% Generate Normal distributions
IndN=find(dist=='N');
if ~isempty(IndN)
    theta(:,IndN)=repmat(mu(IndN)',Nsim,1) + repmat(sigma(IndN)',Nsim,1).*RV(:,IndN);
end
%% Generate Truncated Normal distributions
IndT=find(dist=='T');
if ~isempty(IndT)
    theta(:,IndT)=TruncatedNormal(mu(IndT),sigma(IndT),bound(IndT,1),bound(IndT,2),Nsim);
end

%% Generate Uniform distributions
% mu=a, sigma=b
IndU=find(dist=='U');
if ~isempty(IndU)
    a=mu;
    b=sigma;
    ni=length(IndU);
    theta(:,IndU)= repmat(b(IndU)'-a(IndU)',Nsim,1) .*rand(Nsim,ni)+repmat(a(IndU)',Nsim,1);
end




