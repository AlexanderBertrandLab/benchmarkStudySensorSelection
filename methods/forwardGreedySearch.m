function [groupSel,objFun] = forwardGreedySearch(R1,R2,nbGroupToSel,groupSelector,K)
% FORWARDGREEDYSEARCH Perform a forward group selection search for the 
%   GEVD-problem of the matrix pencil (R1,R2).
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
groupSel = [];
groupRem = 1:nbGroups;
for grpSel = 1:nbGroupToSel
   objFuns = zeros(length(groupRem),1);
   for gr = 1:length(groupRem)
       sel = sum(groupSelector(:,[groupSel,groupRem(gr)]),2);
       E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
       E = sort(E,'descend');
       objFuns(gr) = sum(E(1:min(end,K)))/K; % note: if less filters than required available, use less filters
   end
   [objFun,optGroup] = max(objFuns);
   groupSel = [groupSel,groupRem(optGroup)];
   groupRem(optGroup) = [];
end

end