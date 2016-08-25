function [theta,lnCL_addD,lnpdf_int,count_Change,smcsettings]=...
    MoveSet_MixW(theta, lnCL, int_ComParam, int_Param,lnpdf_int,Y,gamma,smcsettings,CountG,count_Change)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MH move has been separateed by Block
%  1st Block    : Common variables
%  other Blocks : Independent by firms
%  MoveSet__MixW
%  1st block regress on 10 theta from other blocks
%  2nd...Nth block regress on dalta 
%  VAR-COR ( mu, beta, sigma ) 3*3 on each firm
%function assumes an equal-weighted param set

Nfirm=size(theta.(smcsettings.Fname{smcsettings.Ind_Idx(1)}),2);
smcsettings.FigN=smcsettings.FigN+1;
Increm_move=  smcsettings.MoveBlocks(CountG+1).moveindex;  % Incremental part: conditional density

N_indPara = length(smcsettings.Ind_Idx);
Xcondition1=ones(smcsettings.Nsim,1);
for i = 1:N_indPara
    Iname=smcsettings.Ind_Idx(i);
    Xcondition1 = [Xcondition1 theta.(smcsettings.Fname{Iname})];
    Int_data(:,i,:) = squeeze( theta.(smcsettings.Fname{Iname}) );
end
Nreg1 = size(Xcondition1,2)-1;
scale_c=1;
scale_t=1;

AcceptRate=zeros(CountG+1,1);
%loop over Block Moves
for iblock=1:CountG+1
    
    ln_int_Prop = lnpdf_int;
    
    if count_Change(iblock)<smcsettings.ChgB
        
    Block_move = smcsettings.MoveBlocks(iblock).moveindex;
    N1 = Block_move(1);
    N2 = Block_move(end);
    
    theta_proposal=theta;
    
    % Compute proposal density parameters
    Nchg_f=length(Block_move);  % Number of firms Change parameters in iblock
    res = zeros( smcsettings.Nsim ,N_indPara,'double');
       
    if iblock==1; %%% 1st block is common parameter

        %   stats(1)
        rInd = randperm(Nreg1,min(10,Nreg1) );
        Xcon1 = Xcondition1(:,[1 1+rInd]);
        X = theta.(smcsettings.Fname{smcsettings.Com_Idx});
        Reg = ols( X , Xcon1 );
        Xhat = Xcon1 * Reg.beta;
        Rsig2 = cov(Reg.resid);
        
        % q(x,y) = q(x); implements Independent Metropolis-Hastings sampling
        % Generate MultiNormal distribution        
        [X_proposal,logprop_orig,logprop_proposal] = SimProposalMN2( Xhat, Rsig2, X,scale_c, ...
            smcsettings.MoveBlocks(iblock).trcaindx,smcsettings.MoveBlocks(iblock).logisticbound_l, ...
            smcsettings.MoveBlocks(iblock).logisticbound_u);
        
        theta_proposal.(smcsettings.Fname{smcsettings.Com_Idx})=X_proposal;
        % density of initialization sampler for common parmameter
        if CountG==1  % N1 firms
        [ ln_int_Prop(:,1)]=loglikelihood_density( theta_proposal.(smcsettings.Fname{smcsettings.Com_Idx}),...
                int_ComParam.m(1), diag(int_ComParam.s(1)),1,0);
        else  % Adding to Ns firms with mixture delta density function
            ln_int_Prop(:,1) = loglikelihood_density(theta_proposal.(smcsettings.Fname{smcsettings.Com_Idx}),...
                int_ComParam.m, int_ComParam.s,1,1,smcsettings.MixW, int_ComParam.tc);
        end   

        [lnCL ]=feval(smcsettings.fun,theta, Y, 1);
        [l_proposal] = feval(smcsettings.fun ,theta_proposal ,Y ,1 );  

    else
        X_proposal = zeros(smcsettings.Nsim, N_indPara, Nchg_f);
        Xcon2 =[ones(smcsettings.Nsim,1) theta.(smcsettings.Fname{smcsettings.Com_Idx})];
        logprop_proposal=0;
        logprop_orig=0;
        for Gi=1:Nchg_f
            Xhat = zeros(smcsettings.Nsim, N_indPara);
            if Block_move(Gi)==1  % beta_1 without regression
                Xhat = repmat(mean(Int_data(:,:,Block_move(Gi)) ),smcsettings.Nsim,1);
                Xcov = cov(Int_data(:,:,Block_move(Gi)));
            else
            for i = 1: N_indPara
                Reg = ols(Int_data(:,i,Block_move(Gi)),Xcon2);
                res(:,i) = Reg.resid;
                Xhat(:,i) = Xcon2*Reg.beta;
            end
            Xcov = cov(res);
            end

            % q(x,y) = q(x); implements Independent Metropolis-Hastings sampling
            % Generate MultiNormal distribution
            [X_proposal(:,:,Gi),logprop_1,logprop_2] = SimProposalMN2( Xhat, Xcov, Int_data(:,:,Block_move(Gi)),scale_t, ...
                smcsettings.MoveBlocks(iblock).trcaindx(Gi),smcsettings.MoveBlocks(iblock).logisticbound_l, ...
                smcsettings.MoveBlocks(iblock).logisticbound_u);
            logprop_orig=logprop_orig +logprop_1;
            logprop_proposal=logprop_proposal  +logprop_2;
            if CountG==1 || Block_move(Gi)>=Increm_move(1) 
                diagF=1;
                sig2=diag(int_Param(Block_move(Gi)).s);
            [ ln_int_Prop(:,Block_move(Gi)+1)] =loglikelihood_density(...
                    X_proposal(:,:,Gi),int_Param(Block_move(Gi)).m', sig2, diagF,1);
            else
                diagF=0;
            [ ln_int_Prop(:,Block_move(Gi)+1)] =loglikelihood_density(...
                    X_proposal(:,:,Gi),int_Param(Block_move(Gi)).m', int_Param(Block_move(Gi)).s, diagF,1);
            end
        
        end
        for i = 1: N_indPara
            Ifs=smcsettings.Ind_Idx(i);
            theta_proposal.(smcsettings.Fname{Ifs})(:,Block_move) = squeeze(X_proposal(:,i,:));
        end
      
        [l_proposal] = feval(smcsettings.fun,theta_proposal,Y,1);
        [ lnCL] = feval(smcsettings.fun,theta,Y,1);

    end                   
     
    logtargetn_proposal= sum(l_proposal,2)*gamma  + sum(ln_int_Prop,2)*(1-gamma)...
        - logprop_proposal;
    
    logtargetn_orig= sum(lnCL,2 )*gamma + sum(lnpdf_int,2)*(1-gamma)...
        - logprop_orig;


    % log acceptance weights for MH
    % compute acceptance weights
    
    j=find(imag(logtargetn_proposal)~=0);
    
    logtargetn_proposal(j)=-inf;
    
    lnalpha=logtargetn_proposal-logtargetn_orig; %symmetric proposal, difference is zero

    logu=log(rand(smcsettings.Nsim,1));
    
    accept=find(logu<lnalpha);
    
    AcceptRate(iblock)=length(accept)/smcsettings.Nsim;
    
    count_Change(iblock)=count_Change(iblock)+AcceptRate(iblock);
    
    if smcsettings.verbose 
        disp(['Acceptance rate in MH step, group ' num2str(iblock-1) ' : ' num2str(AcceptRate(iblock))]);
    if smcsettings.verbose >1
    figure(mod(smcsettings.FigN,70)+1)
    
%     title(['\gamma= ' num2str(gamma)])
    if iblock==1
    avr_sig1=theta.beta(:,2);
    subplot(1,CountG+1,iblock)
    [sort_S indS]=sort(avr_sig1);
%     scatter(exp(sort_S),theta.weight(indS),'b','.','LineWidth',2)
    scatter(sort_S,theta.weight(indS),'b','.','LineWidth',2)
    title(['weight= ' num2str(mean(theta.weight))])
    xlabel('\beta')
    ylabel('weight')
    hold on
    else
    avr_sig2=theta.beta(:,2);
    subplot(1,CountG+1,iblock)
    [sort_W indW]=sort(theta.weight);
%     scatter(sort_W,exp(avr_sig2(indW)),'b','.','LineWidth',2)
    scatter(sort_W,avr_sig2(indW),'b','.','LineWidth',2)
%     title(['\beta= ' num2str(exp(mean(avr_sig2)))])
    title(['\beta= ' num2str(mean(avr_sig2))])
    ylabel('\beta')
    xlabel('weight')
    hold on
    end
    end
    end
    
    %implement move
    lnpdf_int(accept,:)  = ln_int_Prop(accept,:);
    lnCL_addD=lnCL;
    lnCL_addD(accept,:)  = l_proposal(accept,:);
    if iblock==1
        for i = smcsettings.Com_Idx
            theta.(smcsettings.Fname{i})(accept)=theta_proposal.(smcsettings.Fname{i})(accept);
        end
    else
        tmpDep=0;
%         Xcondition1=ones(smcsettings.Nsim,1);
        for i = smcsettings.Ind_Idx
            theta.(smcsettings.Fname{i})(accept,Block_move)=theta_proposal.(smcsettings.Fname{i})(accept,Block_move);
%             Xcondition1 = [Xcondition1 theta.(smcsettings.Fname{i})];
            tmpDep=tmpDep+Nfirm;
        end
        Int_data(accept,:,Block_move) = X_proposal(accept,:,:);
    end
    
    
    if smcsettings.verbose>1
        if iblock==1
            
            subplot(1,CountG+1,iblock)
%             scatter(exp(sort_S),theta.weight(indS),'g','.','LineWidth',2 )
            scatter(sort_S,theta.weight(indS),'g','.','LineWidth',2 )
            hold off
        else

            avr_sig2=theta.beta(:,2);
            subplot(1,CountG+1,iblock)
%             scatter(sort_W,exp(avr_sig2(indW)) ,'g','.','LineWidth',2 )
            scatter(sort_W,avr_sig2(indW) ,'g','.','LineWidth',2 )
            hold off
        end
    end
    end
    
end
