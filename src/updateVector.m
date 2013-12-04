function update = updateVector(new, current, targetSize)
% update indicator vector

if isempty(current)
    update(1) = new;
else
    n = length(current);
    update = ones(min(n+1,targetSize),1);    
    if n < targetSize
        update(1:n) = current;
        update(n+1) = new;
    else
        for i = 1:targetSize-1;
            update(i) = current(i+1);
        end;
        update(targetSize) = new;
    end;
end;