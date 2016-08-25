% The log-likelihood function of  

% Purpose : Compute log-likelihood function value for a given parameter
% vector

% Input     : PARAM = [mu, sig, beta, weight]  annualized
% data      : n by 7 matrix with the following variables: 
%    1. market capitalization 2. short term debt 3. long term debt 
%    4. other liability 5. total asset 6. interest rate 
%    7. number of trading days between two observatioN1               
% evalue_ini: initial input of the equity value
%          (only for valid data)
% N1: begin of New parameter

function [objval] = Transf_lnL_woMissing(PARAM, data, N1)
% tic
[Nts, Ns]=size(data.Equity);
Npart=size(PARAM.mu,1);

tdays = 250;
   
% preparing parameters for Implied Asset Value
sigma2 = exp( 2*PARAM.nu ) + PARAM.beta.^2;
sig = sqrt( sigma2 );  % Npart * Ns (annumlize)
weight =  PARAM.weight;  

objval = zeros(Npart,1);
sigV = sig(:,1:Ns);  % Note: our model setting Ch_1<=F_n+1

% preparing the log-likelihood function value (only considering 1 year)
diagS2 = exp( 2*PARAM.nu )/tdays;
mu = ( PARAM.mu - sigma2 /2 )/tdays; % Npart * Ns
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


parfor i=1:Npart
    
    % Compute the implied Asset Value
    liability = data.ShortD + 0.5 * data.LongD +  weight(i) * data.otherL; % tdays*Ns
    % To avoid the case that (short term debt + long term debt = 0)
    liability(liability < eps) = eps;
    
    % standardize the Firm Valueby by a Factor of market cap or liability 
    % purpose: fix tolerence level
    Factor = max(data.Equity, liability); % fix at time point 1(1,beg_i:end_i)
    sdEquity = bsxfun(@rdivide,data.Equity,Factor);
    sdStrike = bsxfun(@rdivide,liability,Factor);
    sigVi=sigV(i,:);
    [sdV,Nd1] = impvalue_annual( sdEquity, sdStrike, data.Rf/(100),sigVi , 1, bsxfun(@rdivide,data.TotalA,Factor) );

    impAV = bsxfun(@times, sdV , Factor) ;
%     Nd11=normcdf((log(impAV./liability)+data.Rf(:,:)/100+repmat(sigVi,Nts,1).^2/2)./repmat(sigVi,Nts,1));
    logNd1 = log(Nd1(2:Nts,:)) ;
    
    logVA = log(impAV./data.TotalA);  %%%This part is daily dataccc
    logVAdif = logVA(2:Nts,:)-logVA(1:Nts-1,:); % (tdays-1)*Ns
    
    if N1 == 1     % compute JOINT distr.
        %   Sherman Morrison formula
        invS11 = diag( InvSig2Ns(i,:) )-(AinvbT_Ns(i,:)'*AinvbT_Ns(i,:))/...
            (1+AinvbT_Ns(i,:)*bTNs(i,:)');
        objterm_MN=sum( lnMNpdf( logVAdif, mu1Ns(i,:) ,invS11 ) );
    else  %  Compute the likelihoond associated with CONDITIONAL on preceeding group segment       
        %   Sherman Morrison formula
        invS11s = diag( InvSig2N1(i,:) )-(AinvbT_N1(i,:)'*AinvbT_N1(i,:))/...
            (1+AinvbT_N1(i,:)*bTN1(i,:)');
        CV2= (diag( diagS2(i,1:Ns) )+ bTNs(i,:)'*bTNs(i,:));
        S21=CV2(N1:Ns,1:N1-1);
        S22=CV2(N1:Ns,N1:Ns);
        Y1=logVAdif(:,1:N1-1);
        Y2=logVAdif(:,N1:Ns);
        MUbar=bsxfun(@plus, mu2N1(i,:)' ,S21*invS11s*bsxfun(@minus,Y1,mu1N1(i,:))');
        invSIG2=  inv( S22-S21*invS11s*S21' );
        objterm_MN = sum( lnMNpdf(Y2,MUbar',invSIG2) );
        
    end
    
    TlogVA =  sum(logVA(2:Nts,N1:Ns),2) ;
% Computing the Multivariate Normal distribution log likelihood function
    objval(i,:) = objterm_MN -  sum(TlogVA,1) - sum(sum(logNd1(:,N1:Ns),2) ,1); % sum of time series 
      
end

%  toc

