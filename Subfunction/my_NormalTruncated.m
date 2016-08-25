function sample = my_NormalTruncated(muProposal, varProposal, isConstrained)
% For 1 dim varialbe only, Propose distribution will depend on original mean(muProposal)
% divide into 3 case 1) close to lower B    2) close to Upper bound           
%                    3) in between
% isConstrained(1): lower B  isConstrained(2): Upper B
if ~isscalar( varProposal)
    error('wrong input variant value in proposal stage')
end
STD_pr = sqrt(varProposal);
K = length(muProposal);
pd = makedist('Normal');
sample = zeros(K,1);

% Case close to lower bound
Ind_L=(muProposal-isConstrained(1)<=STD_pr);
if sum( Ind_L )>=1
    MuL= mean(muProposal(Ind_L)); % using average mean
    Lower = (isConstrained(1)-MuL)/ STD_pr;
    Upper = (isConstrained(2)-MuL)/ STD_pr;
    TN = truncate(pd,Lower,Upper);
    sample(Ind_L,1) = MuL + STD_pr*random(TN,sum(Ind_L),1);
end

% Case close to upper bound
SubS = find(~Ind_L);
if ~isempty( SubS )
    Ind_U=(isConstrained(2)-muProposal(SubS)<=STD_pr);
    if sum( Ind_U )>=1
        MuU= mean(muProposal(SubS(Ind_U))); % using average mean
        Lower = (isConstrained(1)-MuU)/ STD_pr;
        Upper = (isConstrained(2)-MuU)/ STD_pr;
        TN = truncate(pd,Lower,Upper);
        sample(SubS(Ind_U),1) = MuU + STD_pr*random(TN,sum(Ind_U),1);
    end
    % Case in between
    SubSubS = SubS(Ind_U~=1);
    while ~isempty(SubSubS)
        sample(SubSubS,1)= my_mvnrnd(muProposal(SubSubS,:),varProposal,1,length(SubSubS));
        SubSubS=find(sample(:,1)>=isConstrained(2) | sample(:,1)<=isConstrained(1));
    end

end