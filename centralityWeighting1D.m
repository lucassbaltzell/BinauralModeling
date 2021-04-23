function [bandWeighted, centralityWeights] = centralityWeighting1D(delayAxis_us, band, CF, weightingFunc, normalization, varargin)
% centralityWeighting - Perform centrality weighting on a 1D array (i.e. specific frequency band).
%
% [BANDWEIGHTED, CENTRALITYWEIGHTS] = centralityWeighting1D(DELAYAXIS_US, CF, BAND, WEIGHTINGFUNC);
%   Perform central weighting on a particular frequency band BAND 
%   with X axis DELAYAXIS_US and characteristic frequency CF, using the
%   function defined in WEIGHTINGFUNC. The weighted band is output in
%   BANDWEIGHTED and its weights in CENTRALITYWEIGHTS.
%
%   Inputs
%
%       DELAYAXIS_US    Array with delay (time) values [us].
%       BAND            Array with values to be weighted.
%       CF              Characteristic frequency [Hz].
%       WEIGHTINGFUNC   String defining the correction method.
%                       Possible values are:
%           'colburn'       Proposed in colburn1977theory (frequency independent).
%           'stern'         Proposed in stern1996lateralization (frequency dependent).
%           'shackleton'    Proposed in shackleton1992across (frequency independent).           
%       NORMALIZATION   String defining what normalization to use.
%                       Possible values are:
%           'pdf'           Values normalized to be a probability density function.
%           'one'           Values normalized so the peak of the weights is one.
%
%   Outputs
%
%       BANDWEIGHTED        Weighted band.
%       CENTRALITYWEIGHTS   Computed and used weights.
%
% Created August 18, 2015.
% Arturo Moncada-Torres


%% Manage optional inputs.
switch numel(varargin)
    case 0
        
    case 1
        addParams = varargin{1};
        if ~mod(length(addParams),2)
            for ii = 1:length(addParams)/2
                paramName   = addParams{(ii*2) - 1};
                paramValue  = addParams{ii*2};
                if ~isempty(paramName) && ~isnan(paramValue)
                    eval([paramName ' = ' num2str(paramValue) ';']);
                end
            end
        end
        
    otherwise
        error([mfilename, ':input'],'Invalid number of additional inputs.');
end


%% Validate inputs.
if length(CF) ~= 1
    error([mfilename,':inputs'],'Invalid number of CF. Only a single CF can be used.');
end
if CF <= 0
    error([mfilename,':inputs'],'Invalid CF value. CF must be larger than zero.');
end


%% Important parameters.
nDelays = length(delayAxis_us);


%% Obtain weights.
switch weightingFunc
    
    case 'stern'
        
        % Define parameters.
        % From stern1996lateralization, Eq. 6 and 7.
        kh              = 3000;
        lf              = 0.1;
        lp              = 1.1;
        if ~exist('flattoplimit_us','var')
            flattoplimit_us = 200;  % [us]
        end

        
        % Eq. 7.
        if CF <= 1200
            kl = lf*(CF^lp);
        else
            kl = lf*(1200^lp);
        end
        
        % Calculate non-normalized weights.
        for delay = 1:nDelays
            
            tau_us = delayAxis_us(delay);   % [us]
            tau_s = tau_us/1e6;             % [us] --> [s]
            
            if abs(tau_us) <= flattoplimit_us
                % The top part of Eq. 6 might be wrong. Instead, use
                % a flat top so that the weight is that at 200 us.
                tau_us = flattoplimit_us;
                tau_s = tau_us/1e6;
            end
            
            exp1 = exp(-2*pi*kl*abs(tau_s));
            exp2 = exp(-2*pi*kh*abs(tau_s));
            weight(delay) = (exp1 - exp2)/abs(tau_s);
        end

        
    case 'colburn'
        
        % Define parameters. Notice how here we use mainly parameters in 
        % ms, in accordance to the original paper.
        % From colburn1977theory, Eq. 3.
        C       = 1;
        limit1	= 0.15;     % [ms]
        limit2	= 2.2;      % [ms]

        % Since this is a frequency-independent weighting function,
        % computation is only done once.
        for delay = 1:nDelays
            
            tau_us = delayAxis_us(delay);   % [us]
            tau_ms = abs(tau_us/1000);      % [ms]
            
            if tau_ms <= limit1
                weight(delay) = C;
                
            elseif (tau_ms > limit1) && (tau_ms <= limit2)
                weight(delay) = C * exp(-(tau_ms - limit1)/0.6);
                
            else
                weight(delay) = C * 0.033 * exp(-(tau_ms - limit2)/2.3);
            end
        end
        
        
    case 'shackleton'
        % Notice that since this is a frequency-independent weighting, it
        % has the same value across all the frequencies.
        
        % Define parameters.
        % From shackleton1992across, p. 2276.
        sd_us = 600;

        % Since this is a frequency-independent weighting function,
        % computation is only done once.
        for delay = 1:nDelays,
            tau_us = delayAxis_us(delay);
            weight(delay) = 1/sqrt((2*pi)*sd_us)  * exp(-1*(tau_us-0)^2/(2*sd_us^2));
        end
        
        
    otherwise
        error('centralityWeighting:inputs','Invalid weightingFunction');     
end

% Perform normalization
switch lower(normalization)
    case 'pdf'
        % PDF normalization (i.e., sum of weights = 1 in each channel).
        C_lf = 1/sum(weight);
        centralityWeights = C_lf .* weight;
    case 'one'
        % One normalization (i.e., make sure max peak of weights is one).
        centralityWeights = weight ./ max(weight);
    otherwise
        error([mfilename,':inputs'],'Invalid normalization method.');
end


%% Apply weights.
bandWeighted = band .* centralityWeights;