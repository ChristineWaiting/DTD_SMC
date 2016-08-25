function [ partData ] = ExtractData( Data, beg_i, end_i, Ncol )
%ExtractData : only output partial Data from beg_i to end_i
% Input structure Data
% row numbers are the same for each field
Fname = fieldnames(Data);
N = length(Fname);
partData = Data;
Erase=[1:beg_i-1 end_i+1:Ncol]; 
for i=1:N
    partData = setfield(partData, Fname{i},{':',Erase},[]);
end
end

