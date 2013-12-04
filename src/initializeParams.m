function params = ...
    initializeParams(indicator, evidenceGathering, stopDecision, params)
% function for checking the params struct and setting default values for
% nonprovided parameters. in case no default value can be set, an error
% message is generated

nPI = length(indicator);
nEGP = length(evidenceGathering);
nSD = length(stopDecision);

% replicate epsilon if a scalar or a PI or EGP-dependent vector is specified
if isfield(params, 'epsilon')
    if length(params.epsilon) == 1
        params.epsilon = repmat(params.epsilon, nEGP, nPI);
    elseif size(params.epsilon, 1) == 1
        params.epsilon = repmat(params.epsilon, nEGP, 1);
    elseif size(params.epsilon, 2) == 1
        params.epsilon = repmat(params.epsilon, 1, nPI);
    elseif ~(size(params.epsilon, 1) == nEGP) && ~(size(params.epsilon, 2) == nPI)
        error('size of the threshold field epsilon must be 1x1, 1xnPI, nEGPx1, or nEGPxnPI');
    end;
else
    % no specification, use recommended values
    epsilon = zeros(nEGP, nPI);
end;

% check input params for the PI and set reasonable threshold values
% you also may set another parameters here (tmem, referenceSet, ...)
for j = 1:nPI
    switch (indicator{j})
        case 'mdr'
            if ~isfield(params, 'epsilon')
                epsilon(:,j) = 0.00002; % Guerrero et al. (CEC 2010)
            end;
        case 'hv'
            if ~isfield(params, 'refPoint')
                error('hypervolumeIndicator function requires a refPoint field in the params struct');
            end;
            if ~isfield(params, 'epsilon')
                epsilon(:,j) = 0.0001; % Wagner et al. (CEC 2010)
            end;
        case 'r'
            if ~isfield(params, 'refPoint') || ~isfield(params, 'idealPoint')
                error('the R indicator requires a refPoint and an idealPoint field in the params struct');
            end;
            if ~isfield(params, 'epsilon')
                epsilon(:,j) = 1e-6; % chosen by Tobias Wagner (no recommendation in a paper)
            end;
        case 'epsilon'
            if ~isfield(params, 'epsilon')
                epsilon(:,j) = 0.0002; % % Guerrero et al. (CEC 2010)
            end;
        case 'cm'
            if ~isfield(params, 'refSet')
                error('convergence metric requires a refSet field in the params struct');
            end;
            if ~isfield(params, 'epsilon')
                epsilon(:,j) = 1e-3; % chosen by Tobias Wagner (no recommendation in a paper)
            end;
        case 'maxCD'
            if ~isfield(params, 'epsilon')
                epsilon(:,j) = 0.02; % Rudenko and Schoenauer (2004)
            end;
            if ~isfield(params, 'tmem')
                params.tmem = 40; % Rudenko and Schoenauer (2004)
            end;
        case 'dqp'
            if ~isfield(params, 'objFunc')
                error('DQP requires a objFunc field in the params struct');
            end;
            if ~isfield(params, 'radius')
                params.radius = 0.05; % Bui et al. (CEC 2009)
            end;
            if ~isfield(params, 'samples')
                params.samples = 500; % Bui et al. (CEC 2009)
            end;
            if ~isfield(params, 'epsilon')
                epsilon(:,j) = 0.01 ; % chosen by Tobias Wagner (no recommendation in a paper)
            end;
        case 'cr'
            % threshold epsilon is set adaptively during the run
            if ~isfield(params, 'tmem')
                params.tmem = 10; % Goel and Stander (IJNME 2010)
            end;
            if ~isfield(params, 'F')
                params.F = 10; % Goel and Stander (IJNME 2010)
            end;
    end;
end;

% check input params for the EGPs and set reasonable default values
for j = 1:nEGP
    switch(evidenceGathering{j})
        case 'direct'
        case 'std'
            if ~isfield(params, 'btrpSamples')
                params.btrpSamples = 10000; % chosen by Tobias Wagner (no recommendation in a paper)
            end;
            if isfield(params, 'tmem') && params.tmem < 2
                error('the standard deviation requires params.tmem > 1');
            end;
        case 'kalman'
            if ~isfield(params, 'R')
                params.R = 0.1; % Marti et al. (CEC 2009)
            end;
        case 'reg'
            if ~isfield(params, 'standardize')
                params.standardize = false;
            end;
            if ~isfield(params, 'combine')
                params.combine = false;
            end;
            if params.combine && ~params.standardize
                warning('standardization is recommended for a combined analysis');
            end;
            if params.standardize && ~isfield(params, 'epsilon')
                epsilon(j,:) = 0; % this is important!
            end;
            if ~isfield(params, 'tmem')
                if params.combine
                    params.tmem = 16; % Wagner et al. (CEC 2010)
                else
                    params.tmem = 30; % Guerrero et al. (CEC 2010)
                end;
            elseif params.tmem < 3
                error('a regression analysis requires params.tmem > 2');
            elseif params.tmem < 10
                warning('params.tmem has a low value. results of regression analysis may be inaccurate');
            end;
        case 'moving'
    end
end

% check input params for the stopping decisions and set reasonable default
% values
for j = 1:nSD
    switch (stopDecision{j})
        case 'threshold'
            if ~isfield(params, 'hits')
                params.hits = 1; % as usual
            end;
        case 'ciNormal'
            if ~isfield(params, 'p')
                params.p = 0.97725; % new recommendation by Luis Marti (no recommendation in a paper)
            end;
            if ~isfield(params, 'hits')
                params.hits = 1; % as usual
            end;
        case 'adaptTest'
            if ~isfield(params, 'p')
                params.p = 0.05; % Wagner et al. (CEC 2010)
            end;
            if ~isfield(params, 'hits')
                params.hits = 2; % Wagner et al. (EMO 2009)
            end;
        case 'validThreshold'
            if ~isfield(params, 'hits')
                params.hits = 1; % Guerrero et al. (CEC 2010)
            end;
    end;
end;

% assign default thresholds
if ~isfield(params, 'epsilon')
    params.epsilon = epsilon;
end;

% default values for necessary parameters if no criterion-based defaults
% were specified
if ~isfield(params, 'tmem')
    params.tmem = 14; % as recommended for OCD-HV by Wagner et al. (CEC 2010)
end;
if ~isfield(params, 'aggregate')
    params.aggregate = 'allPIanyEGP'; % Wagner et al. (EMO 2009)
    % other possibilities are 'any', 'anyPIallEGP', 'all', and 'majority'
end;
if ~isfield(params, 'plot')
    params.plot = true;
end;