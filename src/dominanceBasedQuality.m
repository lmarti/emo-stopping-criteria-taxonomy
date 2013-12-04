function state = dominanceBasedQuality(PFt, state, params)
% implementation of the dominance based quality of set P for the use in 
% the online stopping criteria framework

% DQP works on both objective and decision space
currentSet = PFt.PS;
currentFront = PFt.PF;

% compute dqp for current generation
isDom = false(params.samples,1);
domRatio = zeros(size(currentSet,1),1);
for j = 1:size(currentSet)
    Xtry = randomPointsInCircle(params.samples, currentSet(j,:), ...
            params.radius);
    for k = 1:params.samples                       
        Ftry = feval(params.objFunc, Xtry(k,:));
        ind = paretofront([currentFront(j,:); Ftry]);
        isDom(k) = ~ind(1);        
    end;
    domRatio(j) = mean(isDom);
end;
dqp = mean(domRatio);

% update vector of stored indicator values
if isfield(state, 'dqp')
    state.dqp = updateVector(dqp, state.dqp, params.tmem);
else
    state.dqp = updateVector(dqp, [], params.tmem);
end;