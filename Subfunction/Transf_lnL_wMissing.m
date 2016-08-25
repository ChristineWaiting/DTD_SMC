% The log-likelihood function of  

% Purpose : Compute log-likelihood function value for a given parameter
% vector

% Input     : PARAM = [mu, sig, beta, weight]  annualized
% data1 at time t+1/  data0 : at time t+1
% n by 7 matrix with the following variables: 
%    1. market capitalization 2. short term debt 3. long term debt 
%    4. other liability 5. total asset 6. interest rate 
%    7. number of trading days between two observations               
% evalue_ini: initial input of the equity value
%          (only for valid data)
% N1: begin of New parameter

function [objval] = Transf_lnL_wMissing(PARAM, data1, data0 , N1)
% tic
[Nts, Ns]=size(data1.Equity);
Npart=size(PARAM.mu,1);
tdays = 250;
   
% preparing parameters for Implied Asset Value
sigma2 = exp( 2*PARAM.nu ) + PARAM.beta.^2;
sig = sqrt( sigma2 );  % Npart * Nfirm (annumlize)
weight =  PARAM.weight;  

objval = zeros(Npart,1);
sigV = sig(:,1:Ns);  % Note: our model setting Ch_1<=F_n+1

% preparing the log-likelihood function value (only considering 1 year)
diagS2 = exp( 2*PARAM.nu(:,1:Ns) )/tdays;
mu = ( PARAM.mu(:,1:Ns) - sigma2(:,1:Ns) /2 )/tdays; % Npart * Nfirm
mu1Ns=mu(:,1:Ns);
InvSig2 = 1./diagS2(:,1:Ns);
InvSig2Ns = InvSig2(:,1:Ns);
bTNs= PARAM.beta(:,1:Ns)/sqrt(tdays);
AinvbT_Ns = InvSig2Ns.*bTNs;
% if N1>1
    mu1N1=mu(:,1:N1-1);
    mu2N1=mu(:,N1: Ns);
    InvSig2N1 = InvSig2(:,1:N1-1);
    bTN1= PARAM.beta(:,1:N1-1)/sqrt(tdays);
    AinvbT_N1 = InvSig2N1.*bTN1;
% end

parfor i=1:Npart
    
    % Compute the implied Asset Value
    liability = [data0.ShortD;data1.ShortD] + 0.5 * [data0.LongD;data1.LongD] +...
        weight(i) * [data0.otherL;data1.otherL]; % tdays*Nfirm
    % To avoid the case that (short term debt + long term debt = 0)
    liability(liability < eps) = eps;
    
    % standardize the Firm Valueby by a Factor of market cap or liability 
    % purpose: fix tolerence level
    Factor = max([data0.Equity;data1.Equity], liability); % fix at time point 1(1,beg_i:end_i)
    sdEquity = bsxfun(@rdivide,[data0.Equity;data1.Equity],Factor);
    sdStrike = bsxfun(@rdivide,liability,Factor);
    sigVi=sigV(i,:);
    tmpInt_AV = bsxfun(@rdivide,[data0.TotalA;data1.TotalA],Factor);
    [sdV,Nd1] = impvalue_annual( sdEquity, sdStrike, repmat(data1.Rf/(100),2,1),sigVi , 1, tmpInt_AV );
    
    impAV = bsxfun(@times, sdV , Factor) ;
%     Nd11=normcdf((log(impAV./liability)+data.Rf(:,:)/100+repmat(sigVi,Nts,1).^2/2)./repmat(sigVi,Nts,1));
    logNd1 = log(Nd1(Nts+1:2*Nts,:)) ;
    
    logVA = log(impAV./[data0.TotalA;data1.TotalA]);  %%%This part is daily dataccc
    logVAdif = logVA(Nts+1:2*Nts,:)-logVA(1:Nts,:); % (tdays-1)*Nfirm
    
    %   Sherman Morrison formula
    if N1 == 1     % compute JOINT distr.
    invS11 = diag( InvSig2Ns(i,:) )-(AinvbT_Ns(i,:)'*AinvbT_Ns(i,:))/...
            (1+AinvbT_Ns(i,:)*bTNs(i,:)');
        objterm_MN=sum( lnMNpdf( logVAdif, mu1Ns(i,:) ,invS11 ) ); 
    else  %% F_n < Nfirm 
    invS11 = diag( InvSig2N1(i,:) )-(AinvbT_N1(i,:)'*AinvbT_N1(i,:))/...
            (1+AinvbT_N1(i,:)*bTN1(i,:)');
    %  Compute the likelihoond associated with CONDITIONAL on preceeding group segment
        CV2= (diag( diagS2(i,1:Ns) )+ bTNs(i,:)'*bTNs(i,:));
        S21=CV2(N1:Ns,1:N1-1);
        S22=CV2(N1:Ns,N1:Ns);
        Y1=logVAdif(:,1:N1-1);
        Y2=logVAdif(:,N1:Ns);
        MUbar=bsxfun(@plus, mu2N1(i,:)' ,S21*invS11*bsxfun(@minus,Y1,mu1N1(i,:))');
        invSIG2=  inv( S22-S21*invS11*S21' );
        objterm_MN = sum( lnMNpdf(Y2,MUbar',invSIG2) );
    end
    
    TlogVA =  sum(logVA(Nts+1:2*Nts,N1:Ns),2) ;
% Computing the Multivariate Normal distribution log likelihood function
    objval(i,:) = objterm_MN -  sum(TlogVA,1) - sum(sum(logNd1(:,N1:Ns),2) ,1); % sum of time series 
      
end

%  toc

