function all_seqs = all_synonym_seq(seq_nt, win)
% all_seqs = all_synonym_seq(seq_nt, win)
%   return all possible combinations of synonymous codons in the region
%   given by [win] (in nt, optional, can contain rows for multi-windows).
%   first seq returned is original seq, this ensures that the default 
%   (for equal scores) is the original one. sorting the following
%   seqs by distance from original (does not work for multiple windows).
%
% Alon Diament, Tuller Lab, October 2016.

if nargin < 2
    flank = cell(1, 2);
else
    if min(size(win)) == 1
        % adjust (single) window to fit codons and sequence
        win = [max(0, 3*round(win(1)/3)), 1+min(3*round(win(2)/3), length(seq_nt))];
        if diff(win) < 2
            all_seqs = {seq_nt};
            return
        end
        flank = {seq_nt(1:win(1)), seq_nt(win(2):end)};  % cut
        seq_nt = seq_nt(win(1)+1:win(2)-1);  % paste
    else
        % multiple windows
        % not resorting seqs by distance from origin
        win_combinations = true;
        if win_combinations  % (truly all synonym seqs)
            win = win';
            win(2, :) = win(2, :) + 1;
            all_seqs = arrayfun(@(x, y) {seq_nt(x:y)}, ...
                                [1; win(:)], [win(:)-1; length(seq_nt)], ...
                                'UniformOutput', false);
            all_seqs(2:2:end) = arrayfun(@(x, y) all_synonym_seq(seq_nt(x:y)), ...
                                         win(1, :), win(2, :)-1, 'UniformOutput', false);
            all_seqs = all_combined_seq(all_seqs{:});
            all_seqs = mat2cell(all_seqs, size(all_seqs, 1), ones(1, size(all_seqs, 2)));
            all_seqs = strcat(all_seqs{:});
        else
            % deal with each window independently
            nW = size(win, 1);
            all_seqs = cell(nW, 1);
            for w = 1:nW
                all_seqs{w} = all_synonym_seq(seq_nt, win(w, :));
                if w == 1
                    continue
                end
                all_seqs{w}(1) = [];  % origin
            end
            all_seqs = cat(1, all_seqs{:});
        end
        return
    end
end

if mod(length(seq_nt), 3)
    error('coding seq not divisible by 3');
end

map = revgeneticcode;
seq_aa = nt2aa(seq_nt, 'AlternativeStartCodons', false);
nA = length(seq_aa);

syn = cell(1, nA);
seq_i = zeros(1, nA);
for a = 1:nA
    if seq_aa(a) == '*'
        syn{a} = map.Stops;
    else
        syn{a} = map.(seq_aa(a));
    end
    [~, seq_i(a)] = ismember(seq_nt((a-1)*3 + 1:3*a), syn{a});
end

nS = cellfun(@length, syn);
ncomb = prod(nS);

% all possible combinations
% (using number representation in multiple bases)
all_comb = cell(ncomb, 1);
prior = zeros(ncomb, 1);
comb = zeros(1, nA);
for p = 0:ncomb-1
    tmp = p;
    for b = 1:nA
        comb(b) = mod(tmp, nS(b)) + 1;
        tmp = floor(tmp / nS(b));
    end
    if any(comb == 0)
        continue
    end
    all_comb{p+1} = comb;
end
all_comb = cat(1, all_comb{:});

all_seqs = arrayfun(@(x, y) syn{x}{y}, repmat(1:nA, ncomb, 1), all_comb, 'UniformOutput', false);
all_seqs = arrayfun(@(x) strcat(all_seqs{x, :}), 1:ncomb, 'UniformOutput', false)';
prior = cellfun(@(x) sum(x~=seq_nt), all_seqs);  % distance from origin in nt
all_seqs = strcat(flank{1}, all_seqs, flank{2});

[~, isort] = sort(prior, 'ascend');  % original seq first, followed by similar
all_seqs = all_seqs(isort);

