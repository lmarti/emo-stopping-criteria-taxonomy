function X = randomPointsInCircle(n, c, r)
d = size(c,2);
rad = rand(n,1).*r;
X = ones(n,d).*repmat(rad,1,d);
phi = zeros(n,d-1);
phi(:,1:d-2) = rand(n,d-2).*pi;
phi(:,d-1) = rand(n,1).*2.*pi;
for j = 1:n
    for k = 1:d
        if k > 1
            for l = 1:k-1
                X(j,k) = X(j,k).*sin(phi(j,l));
            end;
        end;
        if k < d
            X(j,k) = X(j,k).*cos(phi(j,k));
        end;
    end
end;
X = repmat(c,n,1) + X;