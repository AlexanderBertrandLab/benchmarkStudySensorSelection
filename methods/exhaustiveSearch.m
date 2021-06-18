function [groupSel,objFun] = exhaustiveSearch(R1,R2,nbGroupToSel,groupSelector,K)
% EXHAUSTIVESEARCH Perform an exhaustive search over all possible group
% selections for the GEVD-problem of the matrix pencil (R1,R2).
%
%   Input parameters:
%       R1 [DOUBLE]: the target covariance matrix
%       R2 [DOUBLE]: the interference covariance matrix
%       nbGroupToSel [INTEGER]: the number of groups to select
%       groupSelector [BINARY]: a nbVariables.nbGroups x nbGroups binary
%           matrix, indicating per group (column) which variables of the 
%           covariance matrices belong to that group with ones at the 
%           corresponding positions.
%       K [INTEGER]: the number of output filters to take into account
%
%   Output parameters:
%       groupSel [INTEGER]: the groups that are selected
%       maxObjFun [DOUBLE]: the corresponding objective (i.e., generalized
%           Rayleigh quotient)

% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

nbGroups = size(groupSelector,2);
cbs = combnk(1:nbGroups,nbGroupToSel); % enumerate all possible combinations
objFuns = zeros(size(cbs,1),1);
for cb = 1:size(cbs,1)
   sel = sum(groupSelector(:,cbs(cb,:)),2);
   E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
   E = sort(E,'descend');
   objFuns(cb) = sum(E(1:min(end,K)))/K; % note: if less filters than required available, use less filters
end
[objFun,optComb] = max(objFuns);
groupSel = cbs(optComb,:)';

end