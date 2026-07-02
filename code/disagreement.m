function d = disagreement(X, N, n)
%DISAGREEMENT  Norma di Frobenius del disaccordo ||x - 1 (x)bar|| nel tempo.
%   X accetta forma (n x N x T)  oppure  (N*n x T).
if ndims(X) == 3
    T = size(X,3); d = zeros(1,T);
    for s = 1:T
        Xi = X(:,:,s); d(s) = norm(Xi - mean(Xi,2), 'fro');
    end
else
    T = size(X,2); d = zeros(1,T);
    for s = 1:T
        Xi = reshape(X(:,s), n, N); d(s) = norm(Xi - mean(Xi,2), 'fro');
    end
end
end
