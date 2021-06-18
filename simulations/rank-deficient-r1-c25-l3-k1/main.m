function main(seed)
rng(str2double(seed))
%% MAIN SCRIPT FOR RAND-DECIFIENT SIMULATIONS
addpath('../..')
addpath('../../methods')
addpath('../../rank-deficiency')

%% Settings and initialization
params = simulationParameters;

if params.nbMCruns > 1
    params.verbose = false;
    params.L1inf.verbose = false;
    params.stecs.verbose = false;
    params.L12.verbose = false;
else
    params.verbose = true;
end

% results structure
results = struct;
results.params = params;
results.methods = params.methods;
chSel = cell(length(params.chRange),length(params.methods),params.nbMCruns);
snr = nan*ones(length(params.chRange),length(params.methods),params.nbMCruns);
compTime = zeros(length(params.chRange),length(params.methods),params.nbMCruns);

% build channel/group selector matrix
params.chSelector = zeros(params.C*params.lags,params.C);
for ch = 1:params.C
    params.chSelector((ch-1)*params.lags+1:ch*params.lags,ch) = 1;
end

%% Monte-Carlo runs
for run = 1:params.nbMCruns
    fprintf('\n Monte-Carlo run %d/%d \n',run,params.nbMCruns);
    
    %% Simulate target and noise signals
    % build individual target signal
    S = 0.5*generateSyntheticData(params.C,params.simulation.Ntarget,params.simulation.T,params.simulation.freqRange,params.simulation.fs,params.simulation.att,params.simulation.maxLag);

    % build time-lagged versions of target signals
    Saug = zeros(params.C*params.lags,params.simulation.T);
    for ch = 1:params.C
        Saug((ch-1)*params.lags+1:ch*params.lags,:) = toeplitz([S(ch,1);zeros(params.lags-1,1)],S(ch,:));
    end
    S = Saug;

    % build noise signal
    N = 6*generateSyntheticData(params.C,params.simulation.Nnoise,params.simulation.T,params.simulation.freqRange,params.simulation.fs,params.simulation.att,params.simulation.maxLag);

    % build time-lagged versions of noise signals
    Naug = zeros(params.C*params.lags,params.simulation.T);
    for ch = 1:params.C
        Naug((ch-1)*params.lags+1:ch*params.lags,:) = toeplitz([N(ch,1);zeros(params.lags-1,1)],N(ch,:));
    end
    N = Naug;

    % build covariance matrices
    Rss = cov(S'); Rnn = cov(N');
    
    %% Test all different methods
    [chSelTemp,snrTemp,compTimeTemp] = testMethods(Rss,Rnn,params);
    chSel(:,:,run) = chSelTemp;
    snr(:,:,run) = snrTemp;
    compTime(:,:,run) = compTimeTemp;
end
results.chSel = chSel;
results.snr = snr;
results.compTime = compTime;

%% Saving results
if params.save
    save('results-'+params.saveName,'results','params');
end

exit
