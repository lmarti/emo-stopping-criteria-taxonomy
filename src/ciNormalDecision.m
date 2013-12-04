function flag = ciNormalDecision(state, params)
% implementation of a threshold decision by means of a normal-distribution
% based confidence interval for the online stopping criteria framework

% initialize flag matrix
flag = zeros(length(state.EGP), length(state.PI));
% initialize memory of former decisions
if ~isfield(state, 'flagPreGen')
    state.flagPreGen = false(length(state.EGP), length(state.PI), ...
        params.hits);
end;

% compute stopping decision
for j = 1:length(state.EGP)
    egpValues = getfield(state, state.EGP{j});
    if ~isnan(egpValues(1,1)) && isnan(egpValues(1,2))
        % no error information for the confidence interval, EGP is ignored
        flag(j,:) = NaN;
        continue;
    end;
    for k = 1:length(state.PI)
        if ~isnan(egpValues(k,1))
            temp = egpValues(k,1) + norminv(params.p).*egpValues(k,2) ...
                < params.epsilon(j,k);
            state.flagPreGen(j,k,:) = updateVector(temp, ...
                reshape(state.flagPreGen(j,k,:), 1, params.hits), ...
                params.hits);
            flag(j,k) = all(state.flagPreGen(j,k,:));
        end;
    end;
end;

% aggregate decision from the set of feasible proposals
flag = flag(~any(isnan(flag),2),:);
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