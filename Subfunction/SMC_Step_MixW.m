function [ESS,normfac,AcceptRate,theta,lnCL_addD,lnw,lnpdf_int,gamma,smcsettings]=...
    SMC_Step_MixW(theta,lnCL_addD,lnw,lnpdf_int, int_ComParam, int_Param, gamma,...
    Y, smcsettings, CountG)
% This SMC_step is for RunSMC_est.m code
%one step of the sequential SMC algorithm
AcceptRate=nan;

normfac=0;

gammavec=linspace(-100,0,1000);

l= sum(lnCL_addD,2) - sum(lnpdf_int,2);

ESSvec=TempFun(gammavec, lnw, l, smcsettings.ESS_bound);

j=find(ESSvec>smcsettings.ESS_bound);

gammaincrement=gammavec(max(j));
    
gammanew=min(gamma+exp(gammaincrement),1);

incrementalw=(gammanew-gamma)*l;

j=find(isnan(incrementalw) | imag(incrementalw)~=0);

incrementalw(j)=-inf;

%get incremental normalizing ratio
W_prev=exp(lnw-max(lnw)); W_prev=W_prev/sum(W_prev);

max_incrementalw=max(incrementalw);

normfac=normfac+log(sum(W_prev.*exp(incrementalw-max_incrementalw)))+max_incrementalw;
%end of normalizing ratio computations

lnw=lnw+incrementalw;

%compute normalized weights and compute ESS
W=exp(lnw-max(lnw));

W=W/sum(W);

ESS=1/sum(W.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if smcsettings.verbose>0
    disp(['Resample at gamma=' num2str(gammanew) ' ESS=' num2str(ESS)]);
end

if smcsettings.verbose >1 && CountG>1
    smcsettings.FigN=smcsettings.FigN+1;
    figure(mod(smcsettings.FigN,70)+1)
    subplot(2,2,1)
    hist(theta.weight)
    title(['\gamma= ' num2str(gamma)])
    xlabel(['weight mean:' num2str(mean(theta.weight))])
    
    avr_sig=mean(theta.beta,2);
    subplot(2,2,3)
    [sort_S, indS]=sort(avr_sig);
    scatter(sort_S,theta.weight(indS),'.')
    %         xlabel(['\sigma mean:' num2str(exp(mean(avr_sig)))])
    xlabel(['\beta mean:' num2str(mean(avr_sig))])
end

[theta,lnw,lnCL_addD,lnpdf_int]=...
    ResampleSet1(theta,lnw,lnCL_addD,lnpdf_int);

if smcsettings.verbose >1 && CountG>1
    subplot(2,2,2)
    hist(theta.weight)
    xlabel(['weight mean:' num2str(mean(theta.weight))])
    
    avr_sig=mean(theta.nu,2);
    subplot(2,2,4)
    [sort_S, indS]=sort(avr_sig);
    scatter(exp(sort_S),theta.weight(indS),'.')
    %         xlabel(['\sigma mean:' num2str(exp(mean(avr_sig)))])
    xlabel(['\nu mean:' num2str(mean(exp(avr_sig)))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Support boosting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ESS<smcsettings.ESS_upper  
       
    counter=0;
    count_Change=zeros(CountG+1,1);
    while ( min(count_Change)<smcsettings.ChgB  ) 
        
        [theta,lnCL_addD,lnpdf_int,count_Change,smcsettings]=...
            MoveSet_MixW(theta,lnCL_addD, int_ComParam, int_Param,lnpdf_int,...
            Y,gammanew,smcsettings, CountG,count_Change);
        
        counter=counter+1;
       
    end

    if smcsettings.verbose
        if gamma==0
            disp('gamma : Acceptant Rate... ')
        end
        disp(['Finish gamma :' num2str(gammanew) ' uisng '  num2str(counter)   ' MH moves ']);
    end
end
gamma=gammanew;








