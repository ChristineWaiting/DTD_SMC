function [y ] = lnMNpdf(X,mu,InvSig2)
% f=-(X-mu).^2./(2.*(sigma.^2))-log(sigma)-log(2*pi)/2;
% X:N*T  % mu: N*T  sigma: N*N
% does not include irrevlent constant
N=size(InvSig2,1);
log_det_Sig2=-( log(det(InvSig2/10000))+ N*log(10000) );
Xm_M=bsxfun(@minus,X,mu);

quadform = sum( Xm_M*InvSig2.*Xm_M ,2 )  ;
logSqrtDetSigma = 0.5*log_det_Sig2 ;
y = bsxfun(@minus,-0.5*quadform , logSqrtDetSigma );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Over 70 random variables may cause det==0 ---> numerical precision issue

% tic
% sum( diag( Xm_M*InvSig2*Xm_M') ) ;
% toc
% Elapsed time is 0.005822 seconds.
% tic
% for i=1:250
% Qi(i)=Xm_M(i,:)*InvSig2*Xm_M(i,:)';
% end
% toc
% Elapsed time is 0.009584 seconds.
% SQ=sum(Qi)
