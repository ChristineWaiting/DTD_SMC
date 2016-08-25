% The log-likelihood function of

% Purpose : Compute log-likelihood function value for a given parameter
% vector with missing data

% Input     : PARAM = [mu, sig, beta, weight]  annualized
% data      : n by 7 matrix with the following variables:
%    1. market capitalization 2. short term debt 3. long term debt
%    4. other liability 5. total asset 6. interest rate
%    7. number of trading days between two observations
% evalue_ini: initial input of the equity value
%          (only for valid data)
% beg_i: begin of firm with parameters changing

function [Outobj ] = Transf_lnL_main(PARAM, data, beg_i)

Fname = fieldnames(data);
nameNF = length(Fname);

Pname = fieldnames(PARAM);
nameNP = length(Pname);

% Missing data treatment
[Nts, Nfirm]=size(data.Equity);
Valid = ~isnan(data.Equity(2:Nts,:)./data.Equity(1:Nts-1,:));
[Gvalid, iE, iG]=unique(Valid,'rows');
N_JD = length(iE);

for i = 1:N_JD
    k=find(iG==i);
    if N_JD > 1
        for Iname = 3:nameNF   %%%% 1:Date 2:FirmCode
            tmp_data1.(Fname{Iname}) = data.(Fname{Iname})(k+1,Gvalid(i,:));
            tmp_data0.(Fname{Iname})=data.(Fname{Iname})(k,Gvalid(i,:));
        end
        tmp_PARAM.(Pname{1}) = PARAM.(Pname{1}); 
        for Iname = 2:nameNP   %%%% 1:Weight 
            tmp_PARAM.(Pname{Iname}) = PARAM.(Pname{Iname})(:,Gvalid(i,:));
        end
        
        [objval(:,i) ] =Transf_lnL_wMissing(tmp_PARAM, tmp_data1,tmp_data0, beg_i );
    else

        [objval(:,i) ] = Transf_lnL_woMissing(PARAM, data, beg_i );
    end

end
Outobj=sum(objval,2) ;