function Wt = maxBlockW(W,groupSelector)
% MAXBLOCKW Compute Wtilde element-wise by taking the maximum over the
% group blocks in W.
%
%   Input parameters:
%       W [DOUBLE]: the extended nbGroups.nbVars x nbGroups.nbVars
%       optimization variable
%       groupSelector [BINARY]: a nbVariables.nbGroups x nbGroups binary
%           matrix, indicating per group (column) which variables of the 
%           covariance matrices belong to that group with ones at the 
%           corresponding positions.

nbGroups = size(groupSelector,2);
Wt = zeros(nbGroups,nbGroups);
for c1 = 1:nbGroups
    sel1 = sum(groupSelector(:,c1),2);
    for c2 = 1:nbGroups
        sel2 = sum(groupSelector(:,c2),2);
        Wt(c1,c2) = max(max(abs(W(sel1==1,sel2==1))));
    end
end