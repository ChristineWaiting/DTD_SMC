function MoveBlocks=DefineBlocksDTD1(N, GroupN,Gnum)
% 1st MHmove Block: Common Variables
MoveBlocks(1).moveindex=1;

MoveBlocks(1).logisticbound_l=0;

MoveBlocks(1).logisticbound_u=1;

MoveBlocks(1).trcaindx=1;

% 2nd MHmove Block: beta>0 (avoid identification issue) 
if Gnum==1
    GroupN = N;
end
MoveBlocks(2).moveindex=1:GroupN;  


MoveBlocks(2).trcaindx=[1 zeros(1,GroupN-1)]; % make sure at least one firm has positve beta

MoveBlocks(2).logisticbound_l=0;

MoveBlocks(2).logisticbound_u=inf;

for i=1:Gnum-1
    
    if i==Gnum-1
        MoveBlocks(i+2).moveindex=( GroupN*i+1:N );
        GN= N-GroupN*i;
    else
        MoveBlocks(i+2).moveindex=( GroupN*i+1 : GroupN*(i+1) );
        GN= GroupN;
    end
    
    MoveBlocks(i+2).logisticbound_l=[];
    
    MoveBlocks(i+2).logisticbound_u=[];
    
    MoveBlocks(i+2).trcaindx=zeros(1,GN);
end
