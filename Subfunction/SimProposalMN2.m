function [X_proposal,lnprop_X,lnprop_Xprop]=SimProposalMN2(Xmean,V_proposal,X, scale,trcaindx,...
    logisticbound_l,logisticbound_u)
%simulate from mv normal fitted on X under constrained upper bound and
%lower bound         input: param_sig2,   output: 
%assume even-weighted sample

[Nparam,K]=size(Xmean);

%fit mv normal on X and sim from it

%transform proposal

X_tr_proposal=my_mvnrnd(Xmean,V_proposal*scale,1,Nparam);

if ( trcaindx == 1 )

    InValid=find(X_tr_proposal(:,1)>=logisticbound_u | X_tr_proposal(:,1)<=logisticbound_l);
    if isscalar( V_proposal )
        X_tr_proposal(InValid,:) = my_NormalTruncated(Xmean(InValid,:), V_proposal*scale, [logisticbound_l logisticbound_u]);
    else
        while ~isempty(InValid)  % keep resample until all beta>0
            %         tmp = my_mvnrnd(X(InValid,:),V_proposal*scale,Nparam,Nparam);
            %         Valid=find(tmp(:,1)<=logisticbound_u | tmp(:,1)>=logisticbound_l);
            %         X_tr_proposal(InValid,:) = tmp( Valid(1:length(InValid)), : );
            X_tr_proposal(InValid,:) = my_mvnrnd(Xmean(InValid,:),V_proposal*scale,1,length(InValid));
            InValid=find(X_tr_proposal(:,1)>=logisticbound_u | X_tr_proposal(:,1)<=logisticbound_l);
        end
    end
end

X_proposal=X_tr_proposal;

% log likelihood of jumping distrubtuion computation
if ~isempty(X)
lnprop_X = loglikelihood_density(X, Xmean',V_proposal ,0 ,1);

lnprop_Xprop = loglikelihood_density(X_proposal, Xmean',V_proposal,0 ,1);
end

end
     


        


