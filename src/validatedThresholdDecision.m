function flag = validatedThresholdDecision(state, params)
% implementation of a threshold decision with an additional evaluation of
% the goodness of fit of the regression analysis for the online stopping
% criteria framework

% initialize flag matrix
flag = false(length(state.EGP), length(state.PI));
% initialize memory of former decisions
if ~isfield(state, 'flagPreGen')
    state.flagPreGen = false(length(state.EGP), length(state.PI), ...
        params.hits);
end;

% compute stopping decision
for j = 1:length(state.EGP)
    egpValues = getfield(state, state.EGP{j});
    for k = 1:length(state.PI)
        if ~isnan(egpValues(k,1))
            if size(egpValues, 2) < 4 || isnan(egpValues(k,4))
                % no residual information -> simple threshold decision
                temp = egpValues(k,1) < params.epsilon(j,k);
            else
                % calculating formula (10) of Guerrero (CEC 2010)
                stabThreshold = 1 - (2./params.tmem) + ...
                    3.*sqrt(2./params.tmem - 4./(params.tmem.^2));
                temp =  (egpValues(k,4) < stabThreshold) && ...
                    (egpValues(k,1) < params.epsilon(j,k));
            end
            state.flagPreGen(j,k,:) = updateVector(temp, ...
                reshape(state.flagPreGen(j,k,:), 1, params.hits), ...
                params.hits);
            flag(j,k) = all(state.flagPreGen(j,k,:));
        end;
    end;
end;

% aggregate decision from the set of proposals
switch params.aggregate
    case 'any'
        flag = any(flag(:));
    case 'allPIanyEGP'
        flag = any(all(flag,2));
    case 'anyPIallEGP'
        flag = all(any(flag,2));
    case 'all'
        flag = all(flag(:));
    case 'majority'
        flag = nansum(flag(:)) > 0.5.*sum(~isnan(flag(:)));
end;