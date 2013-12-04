function state = kalmanEvidenceGathering(state, params )
%   Uses a simplified Kalman filter to determine when to stop.
%
%   Paremeters:
%   * state, the state of the filter in the current iteration. Might
%   be necesary to initialize it at the begining. By default
%   currentState.ind = 1; currentState.P = R.
%   * R, value of the noise of the filter.
%   * indicatorThreshold, a threshold indicating when we should fire.
%
%   Returns:
%   * state, the state of the filter after computation.
%
%   See:
%   Marti, L., Garcia, J., Berlanga, A., & Molina, J. M. (2009).
%   An Approach to Stopping Criteria for Multi-Objective Optimization
%   Evolutionary Algorithms: The MGBM Criterion. In 2009 IEEE Conference
%   on Evolutionary Computation (CEC 2009), Piscataway, New Jersey.
%
%
% ---
% Version: $Id: kalmanEvidenceGathering.m 43 2010-09-17 14:20:14Z wagner $

for j = 1:length(state.PI)
    M = state.(state.PI{j});
    if ~isfield(state, 'kalman') || ...
            size(state.kalman,1) < length(state.PI) || ...
            any(isnan(state.kalman(j,:)))
        % first generation: initialize Kalman filter with first value and
        % initial guess of the noise
        if ~isnan(M(end))
            state.kalman(j,:) = [M(end) params.R 1];
        else
            state.kalman(j,:) = [NaN NaN NaN];
        end;
    elseif ~isnan(M(end))
        % configure the Kalman filter
        s.x = state.kalman(j,1);
        s.z = M(end);
        s.A = 1;
        s.B = 0;
        s.Q = 0;
        s.R = params.R;
        s.P = state.kalman(j,2);
        
        % compute Kalman filter
        out = kalmanf(s);
        
        % update state
        state.kalman(j,:) = [out.x out.P state.kalman(j,3)+1];
    end;
end;