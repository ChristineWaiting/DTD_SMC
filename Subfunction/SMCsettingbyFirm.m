function [ ModelSetting,int_ComParam,int_Param,smcsettings ] = SMCsettingbyFirm( NFirm )
%SMCSETTING : Initial parameters setting for SMC

ModelSetting.N= NFirm; % total number of firm
if ModelSetting.N<3
    error('Error: the minimum company number for computation is 3. ')
end
%  prevent from under identification : N >= 3
ModelSetting.GroupN=5;

%%%% setting parameters
ModelSetting.adjK=1;
smcsettings.Com_Idx=1;
smcsettings.Ind_Idx=[2 3 4];
smcsettings.Fname={'weight','beta','mu','nu'}';
smcsettings.nameN=length(smcsettings.Fname);
ModelSetting.Nc=length(smcsettings.Com_Idx); 
ModelSetting.Ni=length(smcsettings.Ind_Idx);

% Group number
smcsettings.Gnum = ceil(ModelSetting.N / ModelSetting.GroupN);


%% Assume Initial PARAM   
% common parameter: weight
int_ComParam.m=0.5;
int_ComParam.s=0.3.^2;
int_ComParam.dist='T';
int_ComParam.bound=[0 1];
int_ComParam.tc=[];

% individual firm variable: [Beta Mu Nu]
int_Param(1).m = [0.15 0.2 log(0.1)];
int_Param(1).s = (1*[0.05 0.2 0.05]).^2;
int_Param(1).dist = ['T';repmat('N', ModelSetting.Ni-1,1)];  % 1st beta is positive
int_Param(1).bound = [0 inf];

for i=2:ModelSetting.N  
    int_Param(i).m = int_Param(1).m;
    int_Param(i).s = int_Param(1).s;
    int_Param(i).dist = repmat('N',ModelSetting.Ni,1); 
    int_Param(i).bound =[];
end

smcsettings.MixW = 0.8;
smcsettings.Nmoves=20;
smcsettings.verbose=1; % display on command window
smcsettings.MoveBlocks=DefineBlocksDTD1(ModelSetting.N, ModelSetting.GroupN,smcsettings.Gnum);
smcsettings.ChgB=1;
smcsettings.Asymmetric=1;  % 0:Random Walk 1:Independent
smcsettings.fun='Transf_lnL_main';
   
smcsettings.FigN=0;
smcsettings.Nsim=1024; % # of simulating the parameters
smcsettings.ESS_bound=smcsettings.Nsim*.5; % B in paper
smcsettings.ESS_upper=smcsettings.Nsim*0.9;
end

