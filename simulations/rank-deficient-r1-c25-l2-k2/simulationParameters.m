%% PARAMETERS FILE FOR THE RANDOM SIMULATIONS
function params = simulationParameters()
% SIMULATIONPARAMETERS Generate params structure with all parameters for
% the random simulations
%
%   Output parameters:
%       params [STRUCT]: parameter structure

%% general settings
cvx_solver mosek;
params.C = 25; % number of channels
params.lags = 2; % number of lags per channel
params.K = 2; % number of output filters to take into account
params.nbMCruns = 1; % number of (randomized) Monte-Carlo runs
params.chRange = 2:24; % range of channels to select, MUST be increasing for efficiency
params.save = true; % save results or not
params.saveName = 'deficient-r'+string(params.nbMCruns)+'-c'+string(params.C)+'-l'+string(params.lags)+'-k' +string(params.K)+'-'+string(char(floor(24*rand(1, 10))+97)); % name to save results with
[s,git_hash_string] = system('git rev-parse HEAD'); % fetch git hash
params.gitHash = git_hash_string;


%% methods to test
params.methods = {}; % methods to evaluate, comment lines below out/in to select methods
params.methods{end+1} = 'exhaustive search';
params.methods{end+1} = 'random search';
params.methods{end+1} = 'forward search';
params.methods{end+1} = 'backward search';
params.methods{end+1} = 'STECS';
params.methods{end+1} = 'L1,inf-norm selection';
%params.methods{end+1} = 'L1,inf-norm selection (Hamza)';
% params.methods{end+1} = 'SC heuristic';
% params.methods{end+1} = 'L1,2-norm selection';

%% simulating signals settings
params.simulation.T = 1e5; % number of samples to simulate time signals with
params.simulation.fs = 20; % sampling frequency
params.simulation.Ntarget = randi(2*params.C); % number of target signals
params.simulation.Nnoise = randi(2*params.C); % number of noise signals
params.simulation.freqRange = [1,9]; % frequency range to randomly select frequency range per signal from
params.simulation.att = 0.005; % maximal attenuation
params.simulation.maxLag = 2; % maximal lag

%% L1,inf-norm GEVD-based channel selection settings
params.L1inf.lambdaI = 10; % lower bound for binary hyperparameter search
params.L1inf.lambdaLB = 1e-5; % lower bound for binary hyperparameter search
params.L1inf.lambdaUB = 100; % upper bound for binary hyperparameter search
params.L1inf.relTol = 0.1; % relative tolerance to select channels
params.L1inf.nbIt = 15; % number of reweighting iterations
params.L1inf.verbose = true;
params.L1inf.binarySearch.maxIt = 20; % stopping criterion for binary search: if exceeded, NaN solution

%% L1,inf-norm GEVD-based channel selection (HAMZA) settings
params.L1infHamza.lambdaI = 50; % lower bound for binary hyperparameter search
params.L1infHamza.lambdaLB = 1e-5; % lower bound for binary hyperparameter search
params.L1infHamza.lambdaUB = 1e4; % upper bound for binary hyperparameter search
params.L1infHamza.relTol = 0.1; % relative tolerance to select channels
params.L1infHamza.nbIt = 15; % number of reweighting iterations
params.L1infHamza.verbose = true;
params.L1infHamza.binarySearch.maxIt = 20; % stopping criterion for binary search: if exceeded, NaN solution

%% random search settings
params.randomSearch.nbMCruns = 1e3; % number of Monte-Carlo runs for the random search

%% STECS settings
params.stecs.lambdaLB = 0; % lower bound for binary hyperparameter search
params.stecs.lambdaUB = 1e8; % upper bound for binary hyperparameter search
params.stecs.eps = 1e-12; % eps-term to avoid division over zero in gradient
% params.stecs.tol = 0.1; % tolerance to select channels
params.stecs.relTol = 0.1; % relative tolerance to select channels
params.stecs.verbose = true;
params.stecs.binarySearch.maxIt = 100; % stopping criterion for binary search: if exceeded, NaN solution

%% L1,2-norm GEVD-based channel selection settings
params.L12.lambdaLB = 0; % lower bound for binary hyperparameter search
params.L12.lambdaUB = 1e5; % upper bound for binary hyperparameter search
params.L12.tol = 1e-1; % tolerance to remove channels, relative to maximum
params.L12.nbIt = 15; % number of reweighting iterations
params.L12.verbose = false;
params.L12.binarySearch.maxIt = 20; % stopping criterion for binary search: if exceeded, NaN solution

end
