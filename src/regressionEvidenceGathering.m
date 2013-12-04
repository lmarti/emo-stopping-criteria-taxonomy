function state = regressionEvidenceGathering(state, params)
% computing the absolute slope of the values in M as outcome of the
% EGP for the online stopping criteria framework
%
% params
% .standardize: boolean that decides whether the PI are standardized (true)
% before the regression analysis or not (false)
% .combine:     boolean that decides whether all PI are analyzed in a
% combined regression analysis (true) or in separate regression analyses
% (false). for a combined analysis the use of standardization is
% recommended

if params.combine
    n = length(state.PI);
    X = 1:params.tmem;
    % standardization, replicate row vector for each PI
    X = repmat((X-mean(X))./std(X),1,n); 
    Y = zeros(1,n.*params.tmem);
    for j = 1:n
        M = state.(state.PI{j});
        if length(M) < params.tmem || any(isnan(M))
            state.reg = NaN(length(state.PI), 4);
            return;
        else
            PI = M';
        end;        
        if params.standardize && std(PI) > 0 
            PI = (PI-mean(PI))./std(PI); 
        else
            % just centering           
            PI = (PI-mean(PI));
        end;
        % generate row vector of PI
        Y((j-1)*params.tmem+1:j*params.tmem) = PI; 
    end;
    % linear regression without intercept
    betaHat = ((X*X')\X)*Y';
    % compute residuals
    epsilon = Y - X*betaHat;
    % mean squared error of regression
    s2 = (epsilon*epsilon')./(n*params.tmem-1);
    state.reg = repmat([abs(betaHat) sqrt(s2.*inv(X*X')) ...
        n*params.tmem-1 s2], length(state.PI), 4);
else
    for j = 1:length(state.PI)
        M = state.(state.PI{j});
        if length(M) < params.tmem || any(isnan(M))
            state.reg(j,:) = [NaN NaN NaN NaN];
        else
            PI = M';
            X = 1:params.tmem;
            % standardization
            X = (X-mean(X))./std(X);
            if params.standardize && std(PI) > 0
                PI = (PI-mean(PI))./std(PI);
            else
                % just centering
                PI = (PI-mean(PI));
            end;
            % linear regression without intercept
            betaHat = ((X*X')\X)*PI';
            % compute residuals
            epsilon = PI - X*betaHat;
            % mean squared error of regression
            s2 = (epsilon*epsilon')./params.tmem;
            state.reg(j,:) = [abs(betaHat) sqrt(s2.*inv(X*X')) params.tmem s2];
        end;
    end;
end;