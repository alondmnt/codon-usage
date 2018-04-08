function vals = rand_conserv_feat(featdist, epsilon, n_samples)
% vals = rand_conserv_feat(featdist, epsilon, n_samples)
%   given a NxM matrix of distances [featdist] between N sites and M sites,
%   return [n_samples] with replacement from the M set, each within [epsilon]
%   distance from the N set.
%
%   see also: randperm_conserv_feat()
%
% Alon Diament / Tuller Lab, April 2017.

min_cand = 1;  % report # sites with less possible permutations
method = 'uniform';  % uniform / weighted (prefer closer sites)

dim = size(featdist);
assert(length(dim) == 2 && all(featdist(:) >= 0), ...
       '[featdist] must be a non-negative matrix');
dim = dim(1);

if strcmp(method, 'weighted')
    featdist(featdist == 0) = min(featdist(featdist > 0));
end

vals = zeros(dim, n_samples);
cnt = 0;
for samp = 1:n_samples
    % init
    vals(:, samp) = (1:dim)';
    % draw
    cnt = 0;
    for p = 1:dim
        cand = find(featdist(p, :) <= epsilon);
        if length(cand) < min_cand
            [~, cand] = sort(featdist(p, :));
            cand = cand(1:min_cand);  % constant closest neighbors
            cnt = cnt + 1;
        end

        switch method
            case 'uniform'
                vals(p, samp) = cand(randi(length(cand)));
            case 'weighted'
                cdf = cumsum(1./featdist(p, cand) / sum(1./featdist(p, cand)));
                vals(p, samp) = cand(find(cdf >= rand, 1));
        end
    end
end
fprintf('rand_conserv_feat: %d sites >epsilon\n', cnt);
