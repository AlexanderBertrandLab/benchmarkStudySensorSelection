function [groupSel,objFun,lambda,intermediateResults] = l1infSelectionMultOut(R1,R2,nbGroupToSel,groupSelector,params,K)
% L1INFSELECTION Perform L1,inf-norm-based channel selection for the 
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
%           tol [DOUBLE]: tolerance to remove channels, relative to maximum
%           nbIt [INTEGER]: number of reweighting iterations
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
%       intermediateResults [STRUCT ARRAY]: intermediate results in
%           position corresponding to #selected channels

% Author: Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences
% Correspondence: simon.geirnaert@esat.kuleuven.be

%% parameter setting
nbGroups = size(groupSelector,2);
n = size(R1,1);
nbVar = n/nbGroups;
if K > nbGroupToSel*nbVar
    K = nbGroupToSel*nbVar; % limit number of output filters to use by maximum usable
end
nk = n*K;

intermediateResults = repmat(struct('groupSel', nan, 'objFun', nan, 'lambda', nan), nbGroups, 1 );

% initialize tolerance
if rank(R1) < size(R1,1)
    [Vt,Et] = eig(R1);
    [~,ii] = sort(diag(Et),'descend');
    Vt = Vt(:,ii);
    Vt = Vt(:,1:rank(R1));
    R1t = Vt'*R1*Vt; R2t = Vt'*R2*Vt;
    
    [V,E] = eig(R2t,R1t);
    V = Vt*V;
elseif rank(R2) < size(R2,1)
    [Vt,Et] = eig(R2);
    [~,ii] = sort(diag(Et),'descend');
    Vt = Vt(:,ii);
    Vt = Vt(:,1:rank(R2));
    R1t = Vt'*R1*Vt; R2t = Vt'*R2*Vt;
    
    [V,E] = eig(R2t,R1t);
    V = Vt*V;
else
    [V,E] = eig(R2,R1);
end
[~,ii] = min(diag(E));
V = V./diag(V'*R1*V)'; % make sure the eigenvectors are correctly scaled
WtGt = maxBlockW(V(:,ii)*V(:,ii)',groupSelector);
eps = 0.1*std(WtGt(:)); % standard rule: epsilon is 10% of standard deviation expected coeffs

% define regularization parameter relative
q = trace(R2*V(:,ii)*V(:,ii)');

% define tolerance
tolerance = params.relTol*min(abs(diag(WtGt)));

%% channel selection
if nbGroups == nbGroupToSel % if number of groups to select is equal to the total number of groups, solution known
    groupSel = 1:nbGroups;
    lambda = 0;
    nbGroupSel = nbGroupToSel;
else
    %% binary search
    nbGroupSel = -1;
    lambdaLB = params.lambdaLB; lambdaUB = params.lambdaUB;
    itOuter = 1;
    lambda = params.lambdaI;
    while nbGroupSel ~= nbGroupToSel && itOuter <= params.maxIt
        if itOuter > 1
            lambda = lambdaLB+(lambdaUB-lambdaLB)/2;
        end
        B = ones(nbGroups,nbGroups);
        Wtold = zeros(nbGroups,nbGroups);
        Wtnew = Inf(nbGroups,nbGroups);
        it = 1;
        
        while it <= params.nbIt && max(diag(abs(Wtold-Wtnew))) > 0.1*eps
            Wtold = Wtnew;
            % optimization problem
            cvx_begin quiet
                variable W(nk,nk) symmetric;
                variable Wt(nbGroups,nbGroups) symmetric;
                minimize(trace(kron(eye(K),R2)*W)+lambda*q*trace(B*Wt));
                subject to
                    for i = 1:K
                        for j = i:K
                            if i == j
                                trace(R1*W((i-1)*nbGroups*nbVar+1:i*nbGroups*nbVar,(j-1)*nbGroups*nbVar+1:j*nbGroups*nbVar)) == 1;
                            else
                                trace(R1*W((i-1)*nbGroups*nbVar+1:i*nbGroups*nbVar,(j-1)*nbGroups*nbVar+1:j*nbGroups*nbVar)) == 0;
                            end
                        end
                    end
                    W == semidefinite(nk);
                    for ik = 1:K
                        for jk = 1:K
                            for il = 1:nbVar
                                for jl = 1:nbVar
                                    Atemp = abs(W((ik-1)*nbGroups*nbVar+il:nbVar:(ik*nbGroups-1)*nbVar+il,(jk-1)*nbGroups*nbVar+jl:nbVar:(jk*nbGroups-1)*nbVar+jl));
                                    Wt(tril(true(size(Wt)))) >= Atemp(tril(true(size(Atemp))));
                                end
                            end
                        end
                    end
            cvx_end
            
            % iterative reweighting
            B = 1./(Wt+eps);
            % bookkeeping
            Wtnew = Wt;            
            it = it+1;
        end
        
        Wt = maxBlockW(W,groupSelector);
        nbGroupSel = nnz(abs(diag(Wt))>tolerance);
        
        % save intermediate result
        if nbGroupSel > nbGroupToSel && isnan(intermediateResults(nbGroupSel).lambda)
            groupSel = find(abs(diag(Wt))>tolerance);
            sel = sum(groupSelector(:,groupSel),2);
            E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
            E = sort(E,'descend');
            objFun = sum(E(1:min(end,K)))/K; % note: if less filters than required available, use less filters;
            intermediateResults(nbGroupSel).groupSel = groupSel;
            intermediateResults(nbGroupSel).objFun = objFun;
            intermediateResults(nbGroupSel).lambda = lambda;
        end
        

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
        itOuter = itOuter + 1;
    end 
    
    groupSel = find(abs(diag(Wt))>tolerance);
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
    objFun = sum(E(1:min(end,K)))/K; % note: if less filters than required available, use less filters;
end

end
