function meanObjFun = randomSearch(R1,R2,nbGroupToSel,groupSelector,MCruns,K)
% RANDOMSEARCH Perform a random search over possible group
%   selections for the GEVD-problem of the matrix pencil (R1,R2).
%
%   Input parameters:
%       R1 [DOUBLE]: the target covariance matrix
%       R2 [DOUBLE]: the interference covariance matrix
%       nbGroupToSel [INTEGER]: the number of groups to select
%       groupSelector [BINARY]: a nbVariables.nbGroups x nbGroups binary
%           matrix, indicating per group (column) which variables of the
%           covariance matrices belong to that group with ones at the
%           corresponding positions.
%       MCruns [INTEGER]: the number of Monte-Carlo runs
%       K [INTEGER]: the number of output filters to take into account
%
%   Output parameters:
%       meanObjFun [DOUBLE]: the mean objective function (i.e., generalized
%           Rayleigh quotient) obtained over all Monte-Carlo runs

% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

nbGroups = size(groupSelector,2);
objFuns = zeros(MCruns,1);
for run = 1:MCruns
    groupSel = randperm(nbGroups,nbGroupToSel);
    sel = sum(groupSelector(:,groupSel),2);
    E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
    E = sort(E,'descend');
    objFuns(run) = sum(E(1:min(end,K)))/K; % note: if less filters than required available, use less filters;
end
meanObjFun = mean(objFuns);

end