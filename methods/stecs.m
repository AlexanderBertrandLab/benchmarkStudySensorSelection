function [groupSel,objFun,lambda] = stecs(R1,R2,nbGroupToSel,groupSelector,params,K)
% STECS Perform STECS channel selection for the
% GEVD-problem of the matrix pencil (R1,R2).
%
%   Input parameters:
%       R1 [DOUBLE]: the target covariance matrix
%       R2 [DOUBLE]: the interference covariance matrix
%       nbGroupToSel [INTEGER]: the number of groups to select
%       groupSelector [BINARY]: a nbVariables.nbGroups x nbGroups binary
%           matrix, indicating per group (column) which variables of the
%           covariance matrices belong to that group with ones at the
%           corresponding positions.
%       params [STRUCT]: parameter variable, with fields:
%           lambdaLB [DOUBLE]: lower bound for the binary hyperparameter
%                               search
%           lambdaUB [DOUBLE]: upper bound for the binary hyperparameter
%                               search
%           tol [DOUBLE]: tolerance to select channels
%           eps [DOUBLE]: eps-term to avoid division over zero in gradient
%           verbose [BOOLEAN]: display information or not
%           maxIt [INTEGER]: maximimal number of iterations before 
%               conclusion no solution is found
%       K [INTEGER]: the number of output filters to take into account
%
%   Output parameters:
%       groupSel [INTEGER]: the groups that are selected
%       maxObjFun [DOUBLE]: the corresponding objective (i.e., generalized
%           Rayleigh quotient)
%       lambda [DOUBLE]: the hyperparameter at which the group selection
%           was obtained

% REFERENCES TO BE ADDED
% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% parameter settings
% same as in reference
options.Method = 'lbfgs';
options.MaxIter = 1000;
options.optTol = 1e-9;
options.Display = 'off';

% initial solution is solution with all channels
[V,E] = eig(R1,R2);
[~,ii] = max(diag(E));
wInit = V(:,ii);
nbGroups = size(groupSelector,2);
wt = zeros(size(groupSelector,2),1);
for gr = 1:size(groupSelector,2)
    wt(gr) = norm(wInit(groupSelector(:,gr)==1),2);
end
tolerance = params.relTol*min(wt);

%% channel selection
if nbGroups == nbGroupToSel % if number of groups to select is equal to the total number of groups, solution known
    groupSel = 1:nbGroups;
    lambda = 0;
    nbGroupSel = nbGroupToSel;
else
    %% binary search
    nbGroupSel = -1;
    it = 1;
    lambdaLB = params.lambdaLB; lambdaUB = params.lambdaUB;
    while nbGroupSel ~= nbGroupToSel && it <= params.maxIt
        lambda = lambdaLB+(lambdaUB-lambdaLB)/2;
        
        % optimization problem
        objFun = @(w)optFunSTECS(w,R1,R2,lambda,groupSelector,params.eps);
        w = minFunc(objFun,wInit,options);
 
        % select groups
        wt = zeros(size(groupSelector,2),1);
        for gr = 1:size(groupSelector,2)
           wt(gr) = norm(w(groupSelector(:,gr)==1),2);
        end
        nbGroupSel = nnz(wt>tolerance);
        if params.verbose
            fprintf('\n \t Iterating: %d channels selected with lambda %.2e \n',nbGroupSel,lambda);
        end
        
        % hyperparameter update
        if nbGroupSel > nbGroupToSel
            lambdaLB = lambda;
        elseif nbGroupSel < nbGroupToSel
            lambdaUB = lambda;
        end
        
        % update iteration counter
        it = it + 1;
    end
    groupSel = find(wt>tolerance);
end

%% compute objective function
if nbGroupSel ~= nbGroupToSel % no convergence
    groupSel = nan;
    objFun = nan;
    lambda = params.lambdaUB;
else
    sel = sum(groupSelector(:,groupSel),2);
    E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
    E = sort(E,'descend');
    objFun = sum(E(1:min(end,K)))/K;
end

end