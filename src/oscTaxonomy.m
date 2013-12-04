function [stopGeneration state] = oscTaxonomy(algorithm, testcase, nObj, ...
    maxGen, run, indicator, evidenceGathering, stopDecision, params)
% Framework for the evaluation of the online stopping criteria on
% pregenerated data. Each run is stored in a zip
%
% call: [stopGeneration state] = oscTaxonomy(algorithm, testcase, maxGen, run,...
%                                   indicator, evidenceGathering, stopDecision,...
%                                    params)
%
% Input
% algorithm:        name of the moead (string)
% testcase:         name of the test function (string)
% maxGen:           maximum number of generations to test (integer)
% run:              run number (integer)
% indicator:        Pi - PIs for input data generation (cell array of string)
%                   !! transformed for minimization (lower PI are better)
% evidenceGathering:Upsilon - methods to aggregate data (cell array of string)
%                   !! transformed for maximization (lower EGP corresponds to
%                   a lower improvement in the PF)
% stopDecision:     Phi - methods to decide on convergence (cell array of string)
% params:           struct with parameters for indicators, e.g.,
%                   reference or ideal points, the evidence gathering,
%                   e.g., the size of the time window tmem, and the
%                   stopping decision, e.g., the threshold epsilon
%                   (refer to the function initializeParams.m for the
%                    required parameters, their exact names and defaults)
%
% Output
% stopGeneration:   generation for stopping the algorithm
% state:            final state of the criterion
%
% ---
% Version: $Id: runExperiment.m 83 2010-10-05 20:06:36Z marti-big $

% initialize parameters
stopGeneration = NaN;
state.PI = indicator;
state.EGP = evidenceGathering;
params = initializeParams(indicator, evidenceGathering, stopDecision, params);
if params.plot
    figure();
    colors = {'r', 'g', 'b', 'k', 'm', 'c', 'y'};
    hold on;
end;

% preprocess data (can be removed when coupled to an MOEA)
try
    % data already unzipped
    if strcmp(algorithm, 'moead')
        file = sprintf('./%s/%s/%s_%s_%d_2.txt', algorithm, testcase, ...
            algorithm, upper(testcase), run);
    elseif strcmp(algorithm, 'nsga2')
        file = sprintf('./%s/%s/%d/%s_%d_2.txt', algorithm, ...
            testcase, run, testcase, run);
    elseif strcmp(algorithm, 'mocma')
        file = sprintf('./%s/%s/%d/%s_%s_%d_2.txt', algorithm, ...
            testcase, run, algorithm, upper(testcase), run);
    end;
    load(file)
catch ME
    warning(ME.identifier, ME.message);
    % extract archive
    file = sprintf('./%s/%s/%s_%s_%d.zip', algorithm, testcase, ...
        algorithm, testcase, run);
    if ~exist(file, 'file')
        error('OSC:fileNotFound', 'file %s not found', file);
    else
        unzip(file);
    end;
end;

% start generational loop
for i = 1:maxGen
    
    % read PF* (put the interface the MOEA here)
    if strcmp(algorithm, 'moead')
        file = sprintf('./%s/%s/%s_%s_%d_%d.txt', algorithm, testcase, ...
            algorithm, upper(testcase), run, i);
    elseif strcmp(algorithm, 'nsga2')
        file = sprintf('./%s/%s/%d/%s_%d_%d.txt', algorithm, ...
            testcase, run, testcase, run, i);
    elseif strcmp(algorithm, 'mocma')
        file = sprintf('./%s/%s/%d/%s_%s_%d_%d.txt', algorithm, ...
            testcase, run, algorithm, upper(testcase), run, i);
    end;
    if ~exist(file, 'file')
        warning('OSC:fileNotFound', 'data file %s not found', file);
        continue;
    else
        in = load(file);
        opt = paretofront(in(:,1:nObj));
        PFt.PF = in(opt, 1:nObj);
        PFt.PS = in(opt, nObj+1:end);
    end
    
    % compute indicators
    for k = 1:length(indicator) % multiple indicators possible (as in OCD, LSSC)
        switch(indicator{k})
            case 'mdr'
                state = mdrIndicator(PFt, state, params);
            case 'hv'
                state = hypervolumeIndicator(PFt, state, params);
            case 'r'
                state = rInd(PFt, state, params);
            case 'epsilon'
                state = epsInd(PFt, state, params);
            case 'cm'
                state = convergenceMetric(PFt, state, params);
                % case 'dm'
                %   % complex to calculate, will be integrated soon
                %   state = diversityMetric(PFt, state, params);
            case 'maxCD'
                state = maximumCrowdingIndicator(PFt, state, params);
            case 'dqp'
                state = dominanceBasedQuality(PFt, state, params);
            case 'cr'
                state = consolidationRatio(PFt, state, params);
            otherwise
                error('unknown progress indicator!');
        end;
    end;
    
    % update previous PF stored in the state (add your PI if required)
    if ~isempty(strmatch('mdr', indicator, 'exact')) || ...
            ~isempty(strmatch('r', indicator, 'exact')) || ...
            ~isempty(strmatch('epsilon', indicator, 'exact'))
        state = updatePreGen(PFt.PF, state, params);
    end;
    
    % update archives stored in the state (add your PI if required)
    if ~isempty(strmatch('cr', indicator, 'exact'))
        state = updateArchive(PFt.PF, state, params);
    end;
    
    % compute EGP
    for k = 1:length(evidenceGathering) % multiple EGPs possible (as in OCD)
        switch(evidenceGathering{k})
            case 'direct'
                state = directTransfer(state, params);
                if params.plot
                    for j = 1:size(state.direct,1)
                        plot(i, state.direct(j,1), 'Marker', '+', ...
                            'Color', colors{rem(j-1,7)+1});
                    end;
                end;
            case 'std'
                state = stdEvidenceGathering(state, params);
                if params.plot
                    for j = 1:size(state.std,1)
                        plot(i, state.std(j,1), 'Marker', 'o', ...
                            'Color', colors{rem(j-1,7)+1});
                    end;
                end;
            case 'kalman'
                state = kalmanEvidenceGathering(state, params);
                if params.plot
                    for j = 1:size(state.kalman,1)
                        plot(i, state.kalman(j,1), 'Marker', '*', ...
                            'Color', colors{rem(j-1,7)+1});
                    end;
                end;
            case 'reg'
                state = regressionEvidenceGathering(state, params);
                if params.plot
                    for j = 1:size(state.reg,1)
                        plot(i, state.reg(j,1), 'Marker', 'x', ...
                            'Color', colors{rem(j-1,7)+1});
                    end;
                end;
            case 'moving'
                state = movingAverageEvidenceGathering(state, params);
                if params.plot
                    for j = 1:size(state.moving,1)
                        plot(i, state.moving(j,1), 'Marker', 's', ...
                            'Color', colors{rem(j-1,7)+1});
                    end;
                end;
            otherwise
                error('unknown evidence gathering process!');
        end
    end
    
    % decide on stopping the algorithm
    flag = false(length(stopDecision),1);
    for k = 1:length(stopDecision) % multiple stop decisions possible (as in OCD)
        switch(stopDecision{k})
            case 'threshold'
                flag(k) = thresholdDecision(state, params);
            case 'ciNormal'
                flag(k) = ciNormalDecision(state, params);
            case 'adaptTest'
                flag(k) = adaptiveTestingDecision(state, params);
            case 'validThreshold'
                flag(k) = validatedThresholdDecision(state, params);
            otherwise
                error('unknown stopping decision!');
        end;
    end;
    
    % stop algorithm in case all considered stopping criteria hold
    if all(flag)
        stopGeneration = i;
        return;
    end
end;
hold off;