function flag = adaptiveTestingDecision(state, params)
% implementation of a decision function based on EGP-adaptive statistical
% tests for the online stopping criteria framework

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
    if ~isnan(egpValues(1,1)) && ...
            (isnan(egpValues(1,2)) || isnan(egpValues(1,3)))
        % no error information or degrees of freedom for the statistical
        % test, EGP is ignored
        flag(j,:) = NaN;
        continue;
    end;
    switch state.EGP{j}
        case 'std'
            % chi2 variance test for analyzing var(PI) < epsilon
            for k = 1:size(egpValues,1)
                if ~isnan(egpValues(k,1))
                    % compute test statistic
                    Chi = (egpValues(k,1).^2.*egpValues(k,3))./...
                        params.epsilon(j,k).^2;
                    % look up p-value from Chi^2 distribution
                    temp = chi2cdf(Chi, egpValues(k,3)) < params.p;
                    state.flagPreGen(j,k,:) = updateVector(temp, ...
                        reshape(state.flagPreGen(j,k,:), 1, params.hits), ...
                        params.hits);
                    flag(j,k) = all(state.flagPreGen(j,k,:));
                end;
            end;
        otherwise
            % ttest test for analyzing mean(PI) > epsilon
            for k = 1:size(egpValues,1)
                if ~isnan(egpValues(k,1))
                    % compute test statistic
                    t = (egpValues(k,1)-params.epsilon(j,k))./ ...
                        (egpValues(k,2)./egpValues(k,3)+1);
                    % look up p-value from t-student distribution
                    temp = 2*min(tcdf(t, egpValues(k,3)), ...
                        1-tcdf(t, egpValues(k,3))) > params.p;
                    state.flagPreGen(j,k,:) = updateVector(temp, ...
                        reshape(state.flagPreGen(j,k,:), 1, params.hits), ...
                        params.hits);
                    flag(j,k) = all(state.flagPreGen(j,k,:));
                end;
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