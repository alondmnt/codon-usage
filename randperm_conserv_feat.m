function perm = randperm_conserv_feat(featdist, epsilon, n_samples)
% perm = randperm_conserv_feat(featdist, epsilon, n_samples)
%   given a matrix of distances [featdist] between sites, return [n_samples] 
%   permuations that preserve these distances within a radius [epsilon].
%
% Alon Diament / Tuller Lab, April 2017.

min_cand = 2;  % report # sites with less possible permutations
method = 'uniform';  % uniform / weighted (prefer closer sites)

dim = size(featdist);
assert(length(dim) == 2 && dim(1) == dim(2) && all(featdist(:) >= 0), ...
       '[featdist] must be a non-negative square matrix');
dim = dim(1);

if strcmp(method, 'weighted')
    featdist(featdist == 0) = min(featdist(featdist > 0));
    featdist(1:dim+1:dim^2) = epsilon;  % include self, non-zero
end

perm = zeros(dim, n_samples);
cnt = 0;
for samp = 1:n_samples
    % init
    perm(:, samp) = (1:dim)';
    % permute
    for i = 1:5
        cnt = 0;
        for p = 1:dim
            cand = find(featdist(p, :) <= epsilon);
            cand = cand(featdist(cand, p) <= epsilon);  % vice versa (important!)
            if length(cand) < min_cand
%                 [~, cand] = sort(featdist(p, :));
%                 cand = cand(1:min_cand);  % constant closest neighbors
                cnt = cnt + 1;
            end

            switch method
                case 'uniform'
                    next = cand(randi(length(cand)));
                case 'weighted'
                    cdf = cumsum(1./featdist(p, cand) / sum(1./featdist(p, cand)));
                    next = cand(find(cdf >= rand, 1));
            end

            perm([p, next], samp) = perm([next, p], samp);

            % this ensures that [featdist] follows perm and next neighbors are
            % selected correctly. permuting only columns to keep distance from
            % *original* site (before permutation).
            featdist(:, [p, next]) = featdist(:, [next, p]);
        end
    end
    featdist(:, perm(:, samp)) = featdist;  % permute back to origin
end
fprintf('randperm_conserv_feat: %d sites not permuted\n', cnt);
