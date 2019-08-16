% Finds Euclidean norm between vectors in two array.
function dist = Distance_Function (V, X)
c = size(V, 1);
p = size(V, 2);
n = size(X, 1);
dist = zeros(c, n);
if p > 1,
    for k = 1:c,
	dist(k, :) = sqrt( sum(((X-ones(n, 1)*V(k, :)).^2)') );
    end
else	% 1-D data
    for k = 1:c,
	dist(k, :) = abs(V(k)-X)';
    end
end
