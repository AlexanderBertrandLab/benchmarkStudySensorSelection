function X = generateSyntheticData(C,N,T,freqRange,fs,att,maxLag)
% GENERATESYNTHETICDATA Generate synthetic (EEG) data with a point-source
% model. Sensors are positioned on a square grid between (1,1) and
% (sqrt(C),sqrt(C)), sources can be located in a grid between (0,0) and 
% (sqrt(C)+1,sqrt(C)+1).
%
%   Input parameters:
%       C [INTEGER]: the number of sensors, square number
%       N [INTEGER]: the number of sources
%       T [INTEGER]: the number of time samples
%       freqRange [DOUBLE 2x1]: the lower and upper frequencies to filter
%       randomly within
%       fs [DOUBLE]: the sampling frequency
%       att [DOUBLE]: the maximal attentuation (relative to 1), smaller than
%           one
%       maxLag [INTEGER]: the maximal lag at the furthest possible distance
%           between sensor and source (in samples)
%
%   Output parameters:
%       X [DOUBLE CxT]: the simulated synthetic data

% Authors: 
% Jonathan Dan, KU Leuven, ESAT & Byteflies
% Simon Geirnaert, KU Leuven, ESAT & Dept. of Neurosciences

%% parameter handling
if sqrt(C) ~= round(sqrt(C))
    error('C is expected to be a square number');
end
spread = sqrt(C/log(1/att)); % compute spread such that maximal attenuation is = p

%% source generation and spreading
X = zeros(C,T);
for sc = 1:N
   location = rand(2,1)*(sqrt(C)+1); % location of source 
   
   % source generation
   hp = randi([freqRange(1),freqRange(2)-1]);
   lp = randi([hp+1,freqRange(2)]);
   d = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',hp,'HalfPowerFrequency2',lp,'SampleRate',fs);
   s = randn(1,T); s = (1+0.3*randn(1,1))*filtfilt(d,s);
   
   % source spreading
   for x = 1:sqrt(C)
       for y = 1:sqrt(C)
            distance = sqrt((location(1)-x)^2 + (location(2)-y)^2);
            lag = round(distance/sqrt(2*C)*maxLag);
            X(x+(y-1)*sqrt(C),:) = X(x+(y-1)*sqrt(C),:) + exp(-distance^2/(2*spread^2)).*circshift(s,lag);
       end
   end
end

%% instrumentation noise
X = X + 2*att*randn(C,T);

end