function [chSelected,snr,compTime] = testMethods(R1,R2,params)
% TESTMETHODS Test different channel selection methods on the GEVD problem
% with the matrix pencil (R1,R2).
%
%   Input parameters:
%       R1 [DOUBLE]: the target covariance matrix
%       R2 [DOUBLE]: the noise covariance matrix
%       params [STRUCT]: the parameter structure with all info
%
%   Output parameters:
%       chSelected [CELL chRange x nbMethods]: the selected channels per 
%           number of channels and method
%       snr [DOUBLE chRange x nbMethods]: the corresponding SNRs
%       compTime [DOUBLE chRange x nbMethods]: the corresponding required
%           computation times

% Authors: 
% Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences

%% Initialization
chSelected = cell(length(params.chRange),length(params.methods));
snr = zeros(length(params.chRange),length(params.methods));
compTime = zeros(length(params.chRange),length(params.methods));

%% Exhaustive search
if any(cellfun(@(meth)strcmp(meth,'exhaustive search'),params.methods))
    if params.verbose
        fprintf('\n%s\n*** Exhaustive search ***\n%s\n',repmat('-',1,30),repmat('-',1,30))
    end
    ind = find(cellfun(@(meth)strcmp(meth,'exhaustive search'),params.methods));
    for chr = 1:length(params.chRange)
        if params.verbose
            fprintf('\n %d channels selecting \n',params.chRange(chr));
        end
        tic;
        [chSel,objFun] = exhaustiveSearch(R1,R2,params.chRange(chr),params.chSelector,params.K);
        compTimeTemp = toc;
        chSelected{chr,ind} = chSel; snr(chr,ind) = 10*log10(objFun); compTime(chr,ind) = compTimeTemp;
    end
end

%% Random search
if any(cellfun(@(meth)strcmp(meth,'random search'),params.methods))
    if params.verbose
        fprintf('\n%s\n*** Random search ***\n%s\n',repmat('-',1,30),repmat('-',1,30))
    end
    ind = find(cellfun(@(meth)strcmp(meth,'random search'),params.methods));
    for chr = 1:length(params.chRange)
        if params.verbose
            fprintf('\n %d channels selecting \n',params.chRange(chr));
        end
        tic;
        meanObjFun = randomSearch(R1,R2,params.chRange(chr),params.chSelector,params.randomSearch.nbMCruns,params.K);
        compTimeTemp = toc;
        snr(chr,ind) = 10*log10(meanObjFun); compTime(chr,ind) = compTimeTemp;
    end
end

%% Forward greedy search
if any(cellfun(@(meth)strcmp(meth,'forward search'),params.methods))
    if params.verbose
        fprintf('\n%s\n*** Forward search ***\n%s\n',repmat('-',1,30),repmat('-',1,30))
    end
    ind = find(cellfun(@(meth)strcmp(meth,'forward search'),params.methods));
    for chr = 1:length(params.chRange)
        if params.verbose
            fprintf('\n %d channels selecting \n',params.chRange(chr));
        end
        tic;
        [chSel,objFun] = forwardGreedySearch(R1,R2,params.chRange(chr),params.chSelector,params.K);
        compTimeTemp = toc;
        chSelected{chr,ind} = chSel; snr(chr,ind) = 10*log10(objFun); compTime(chr,ind) = compTimeTemp;
    end
end

%% Backward greedy search
if any(cellfun(@(meth)strcmp(meth,'backward search'),params.methods))
    if params.verbose
        fprintf('\n%s\n*** Backward search ***\n%s\n',repmat('-',1,30),repmat('-',1,30))
    end
    ind = find(cellfun(@(meth)strcmp(meth,'backward search'),params.methods));
    for chr = 1:length(params.chRange)
        if params.verbose
            fprintf('\n %d channels selecting \n',params.chRange(chr));
        end
        tic;
        [chSel,objFun] = backwardGreedySearch(R1,R2,params.chRange(chr),params.chSelector,params.K);
        compTimeTemp = toc;
        chSelected{chr,ind} = chSel; snr(chr,ind) = 10*log10(objFun); compTime(chr,ind) = compTimeTemp;
    end
end

%% STECS
if any(cellfun(@(meth)strcmp(meth,'STECS'),params.methods))
    if params.verbose
        fprintf('\n%s\n*** STECS ***\n%s\n',repmat('-',1,45),repmat('-',1,45))
    end
    ind = find(cellfun(@(meth)strcmp(meth,'STECS'),params.methods));
    
    paramsSTECS = struct;
    paramsSTECS.lambdaLB = params.stecs.lambdaLB;
    lambdaPrev = params.stecs.lambdaUB;
    paramsSTECS.relTol = params.stecs.relTol;
    paramsSTECS.eps = params.stecs.eps;
    paramsSTECS.verbose = params.stecs.verbose;
    paramsSTECS.maxIt = params.stecs.binarySearch.maxIt;
    
    for chr = 1:length(params.chRange)
        if params.verbose
            fprintf('\n %d channels selecting \n',params.chRange(chr));
        end
        paramsSTECS.lambdaUB = lambdaPrev;
        tic;
        [chSel,objFun,lambdaPrev] = stecs(R1,R2,params.chRange(chr),params.chSelector,paramsSTECS,params.K);
        compTimeTemp = toc;
        % add random previous channel if no solution is found
        if any(isnan(chSel)) && chr > 1 && ~any(isnan(chSelected{chr-1,ind}))
            if params.verbose
                fprintf('\n No convergence, resorting to n-1 channels + random \n');
            end
            chSel = chSelected{chr-1,ind};
            notSelected = setdiff(1:length(params.chRange),chSel);
            chSel = sort([chSel;notSelected(randi(length(notSelected)))],'ascend'); % select random channel
            sel = sum(params.chSelector(:,chSel),2);
            E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
            E = sort(E,'descend');
            objFun = sum(E(1:min(end,params.K)));
        end
        chSelected{chr,ind} = chSel; snr(chr,ind) = 10*log10(objFun); compTime(chr,ind) = compTimeTemp;
    end
    
end

%% L1,inf-norm GEVD-based selection
if any(cellfun(@(meth)strcmp(meth,'L1,inf-norm selection'),params.methods))
    if params.verbose
        fprintf('\n%s\n*** L1,inf-norm GEVD-based selection ***\n%s\n',repmat('-',1,45),repmat('-',1,45))
    end
    ind = find(cellfun(@(meth)strcmp(meth,'L1,inf-norm selection'),params.methods));
    
    paramsL1inf = struct;
    paramsL1inf.lambdaLB = params.L1inf.lambdaLB;
    lambdaPrev = params.L1inf.lambdaUB;
    paramsL1inf.lambdaI = params.L1inf.lambdaI;
    paramsL1inf.relTol = params.L1inf.relTol;
    paramsL1inf.nbIt = params.L1inf.nbIt;
    paramsL1inf.verbose = params.L1inf.verbose;
    paramsL1inf.maxIt = params.L1inf.binarySearch.maxIt;
    intermediateLambda = zeros(length(params.chRange),1);
    
    for chr = 1:length(params.chRange)
        if params.verbose
            fprintf('\n %d channels selecting \n',params.chRange(chr));
        end
        if isempty(chSelected{chr,ind})
            paramsL1inf.lambdaUB = lambdaPrev;
            if chr > 1
            	paramsL1inf.lambdaI = paramsL1inf.lambdaLB+(paramsL1inf.lambdaUB-paramsL1inf.lambdaLB)/2;
            end
            if any(intermediateLambda(chr+1:end)~=0)
                nz = nonzeros(intermediateLambda(chr+1:end));
                paramsL1inf.lambdaLB = nz(1);
            else
                paramsL1inf.lambdaLB = params.L1inf.lambdaLB;
            end
            tic;
            [chSel,objFun,lambdaPrev,intermediateResults] = l1infSelectionMultOut(R1,R2,params.chRange(chr),params.chSelector,paramsL1inf,params.K);
            compTimeTemp = toc;
            % add random previous channel is no solution is found
            if any(isnan(chSel)) && chr > 1 && ~any(isnan(chSelected{chr-1,ind}))
                if params.verbose
                    fprintf('\n No convergence, resorting to n-1 channels + random \n');
                end
                chSel = chSelected{chr-1,ind};
                notSelected = setdiff(1:length(params.chRange),chSel);
                chSel = sort([chSel;notSelected(randi(length(notSelected)))],'ascend'); % select random channel
                sel = sum(params.chSelector(:,chSel),2);
                E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
                E = sort(E,'descend');
                objFun = sum(E(1:min(end,params.K)));
            end
            chSelected{chr,ind} = chSel; snr(chr,ind) = 10*log10(objFun); compTime(chr,ind) = compTimeTemp;
            % store intermediate results
            for chrIn = chr+1:length(params.chRange)
                if ~isnan(intermediateResults(params.chRange(chrIn)).groupSel)
                    chSelected{chrIn,ind} = intermediateResults(params.chRange(chrIn)).groupSel;
                    snr(chrIn,ind) = 10*log10(intermediateResults(params.chRange(chrIn)).objFun);
                    compTime(chrIn,ind) = 0;
                    intermediateLambda(chrIn) = intermediateResults(params.chRange(chrIn)).lambda;
                end
            end
        end
    end
end

%% L1,inf-norm GEVD-based selection (Hamza version)
if any(cellfun(@(meth)strcmp(meth,'L1,inf-norm selection (Hamza)'),params.methods)) && params.K > 1
    warning('STECS cannot be combined with more than one filter ouput');
elseif any(cellfun(@(meth)strcmp(meth,'L1,inf-norm selection (Hamza)'),params.methods))
    if params.verbose
        fprintf('\n%s\n*** L1,inf-norm GEVD-based selection (Hamza) ***\n%s\n',repmat('-',1,45),repmat('-',1,45))
    end
    ind = find(cellfun(@(meth)strcmp(meth,'L1,inf-norm selection (Hamza)'),params.methods));
    
    paramsL1infHamza = struct;
    paramsL1infHamza.lambdaLB = params.L1infHamza.lambdaLB;
    paramsL1infHamza.lambdaI = params.L1infHamza.lambdaI;
    lambdaPrev = params.L1infHamza.lambdaUB;
    paramsL1infHamza.relTol = params.L1infHamza.relTol;
    paramsL1infHamza.nbIt = params.L1infHamza.nbIt;
    paramsL1infHamza.verbose = params.L1infHamza.verbose;
    paramsL1infHamza.maxIt = params.L1infHamza.binarySearch.maxIt;
    intermediateLambda = zeros(length(params.chRange),1);
    
    for chr = 1:length(params.chRange)
        if params.verbose
            fprintf('\n %d channels selecting \n',params.chRange(chr));
        end
        if isempty(chSelected{chr,ind})
            paramsL1infHamza.lambdaUB = lambdaPrev;
            if chr > 1
                paramsL1infHamza.lambdaI = paramsL1infHamza.lambdaLB+(paramsL1infHamza.lambdaUB-paramsL1infHamza.lambdaLB)/2;
            end
            if any(intermediateLambda(chr+1:end)~=0)
                nz = nonzeros(intermediateLambda(chr+1:end));
                paramsL1infHamza.lambdaLB = nz(1);
            else
                paramsL1infHamza.lambdaLB = params.L1infHamza.lambdaLB;
            end
            tic;
            [chSel,objFun,lambdaPrev,intermediateResults] = l1infSelectionMultOutHamza(R1,R2,params.chRange(chr),params.chSelector,paramsL1infHamza);
            compTimeTemp = toc;
            % add random previous channel is no solution is found
            if any(isnan(chSel)) && chr > 1 && ~any(isnan(chSelected{chr-1,ind}))
                if params.verbose
                    fprintf('\n No convergence, resorting to n-1 channels + random \n');
                end
                chSel = chSelected{chr-1,ind};
                notSelected = setdiff(1:length(params.chRange),chSel);
                chSel = sort([chSel;notSelected(randi(length(notSelected)))],'ascend'); % select random channel
                sel = sum(params.chSelector(:,chSel),2);
                E = eig(R1(sel==1,sel==1),R2(sel==1,sel==1));
                E = sort(E,'descend');
                objFun = sum(E(1:min(end,params.K)));
            end
            chSelected{chr,ind} = chSel; snr(chr,ind) = 10*log10(objFun); compTime(chr,ind) = compTimeTemp;
            % store intermediate results
            for chrIn = chr+1:length(params.chRange)
                if ~isnan(intermediateResults(params.chRange(chrIn)).groupSel)
                    chSelected{chrIn,ind} = intermediateResults(params.chRange(chrIn)).groupSel;
                    snr(chrIn,ind) = 10*log10(intermediateResults(params.chRange(chrIn)).objFun);
                    compTime(chrIn,ind) = 0;
                    intermediateLambda(chrIn) = intermediateResults(params.chRange(chrIn)).lambda;
                end
            end
        end
    end
end

end