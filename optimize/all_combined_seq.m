function all_seqs = all_combined_seq(varargin)
% all_seqs = all_combined_seq(varargin)
%   return all possible combinations of sequences, one from each of the
%   options listed in each of the arguments. the resulting sequences will
%   be a concatenation of [nargin] parts.
%   returned seqs are sorted by distance from original sequence.
%
% Alon Diament, Tuller Lab, October 2016.

nS = cellfun(@length, varargin);
L = length(varargin);
ncomb = prod(nS);

% all possible combinations
% (using number representation in multiple bases)
all_comb = cell(ncomb, 1);
comb = zeros(1, L);
for p = 0:ncomb-1
    tmp = p;
    for b = 1:L
        comb(b) = mod(tmp, nS(b)) + 1;
        tmp = floor(tmp / nS(b));
    end
    if any(comb == 0)
        continue
    end
    all_comb{p+1} = comb;
end
all_comb = cat(1, all_comb{:});

all_seqs = arrayfun(@(x, y) varargin{x}{y}, repmat(1:L, ncomb, 1), all_comb, 'UniformOutput', false);

% slower, accurate sorting
prior = strcat(all_seqs(:, 1), all_seqs(:, 2));
prior = cellfun(@(x) sum(x~=prior{1}), prior);  % distance from origin in nt

[~, isort] = sort(prior, 'ascend');  % original seq first, followed by similar
all_seqs = all_seqs(isort, :);

