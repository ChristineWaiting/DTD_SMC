function [ lnpdf, trunConst ]=loglikelihood_density(theta, mu, varcov, diagInd, NTind, MixW,trunConst)
%
% diagInd =1: varcov is diag matrix
%         =0: varcov is non-diag matrix
% truncated normal has same densitry as normal except for irrivelent canstant.
% NTind =1 : normal distr.   (for initial sampler)
%       =0 : truncated normal  (for prior sampler)
[Nparam,k]=size(theta);
lnpdf=zeros(Nparam,k);
if nargin==5
    
        if varcov==0
            lnpdf = zeros(Nparam,1);
        else
            MUbar   = mu';
            if diagInd ==1
                Dsig = diag(varcov);
                invSIG2 = diag( 1./Dsig );
            else
                invSIG2 = inv(varcov);
            end
            lnpdf= lnMNpdf(theta,MUbar,invSIG2 );
        end
    
    if NTind==0 % truncated Normal
        %   lnpdf(:,Uind)=repmat(-log([prior_param(Uind).sigma]-[prior_param(Uind).mu]),Nparam,1);
        sigma=varcov.^(1/2);
        trunConst=log(sigma*(normcdf(1,mu',sigma)-normcdf(0,mu',sigma)));
        lnpdf=lnpdf-trunConst;
    end
else % Mixed Truncated Normal
    % mu1: initial setting (1-MixW)
    % mu2: random by last posterior (MixW)
    sigma=varcov.^(1/2);
    if isempty(trunConst)
    % Need constant as Common Variables are constrained to [0,1]
    trunConst(1)=sigma(1)*(normcdf(1,mu(1),sigma(1))-normcdf(0,mu(1),sigma(1))); %I_0
    trunConst(2)=sigma(2)*(normcdf(1,mu(2),sigma(2))-normcdf(0,mu(2),sigma(2))); %I_p
    end
    LnP1=lnMNpdf(theta,mu(1),sigma(1)*sigma(1));
    LnP2=lnMNpdf(theta,mu(2),sigma(2)*sigma(2));
    Const12=repmat(max(max(LnP1,LnP2)),Nparam,1) ;
    lnpdf=log( (1-MixW)*exp(LnP1-Const12)/trunConst(1)+MixW*exp(LnP2-Const12)/trunConst(2) )+Const12;
    
end


