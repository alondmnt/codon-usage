function all_seqs = all_nt_seq(seq_nt, win, max_comb)
% all_seqs = all_nt_seq(seq_nt, win)
%   return all possible combinations of synonymous codons in the region
%   given by [win] (in nt, optional, can contain rows for multi-windows).
%   first seq returned is the original seq, this ensures that the default 
%   (for equal scores) is the original one. sorting the following
%   seqs by distance from original (not implemented for multiple windows).
%
% Alon Diament, Tuller Lab, Oct 2016.

if nargin < 3
    max_comb = 4^9;  % when # seqs exceeds, split sequence (lose combinations)
end
if nargin < 2
    flank = cell(1, 2);
else
    if min(size(win)) == 1
        % adjust (single) window to fit sequence
        win = [max(1, win(1))-1, 1+min(win(2), length(seq_nt))];
        if diff(win) < 2
            all_seqs = {seq_nt};
            return
        end
        flank = {seq_nt(1:win(1)), seq_nt(win(2):end)};  % cut
        seq_nt = seq_nt(win(1)+1:win(2)-1);  % paste
    else
        % multiple windows
        % deal with each window independently
        % not resorting seqs by distance from origin
        nW = size(win, 1);
        all_seqs = cell(nW, 1);
        for w = 1:nW
            all_seqs{w} = all_nt_seq(seq_nt, win(w, :));
            if w == 1
                continue
            end
            all_seqs{w}(1) = [];  % origin
        end
        all_seqs = cat(1, all_seqs{:});
        return
    end
end

subs = {'A', 'C', 'G', 'T'};
nS = length(subs);
L = length(seq_nt);
ncomb = nS^L;

if ncomb > max_comb
    % protection mechanism to handle long windows
    fprintf('all_nt_seq: split\n');
    mid = round(L/2);
    all_seqs = all_nt_seq(seq_nt, [1 mid], max_comb / 2);
    erase = size(all_seqs, 1) + 1; % first seq returned is original seq
    all_seqs = [all_seqs;
                all_nt_seq(seq_nt, [mid+1, L], max_comb / 2)];
    all_seqs = strcat(flank{1}, all_seqs([1:erase-1, erase+1:end]), flank{2});
    return
end

seq_i = nt2int(seq_nt);

% all possible combinations
% (using number representation)
all_comb = cell(ncomb, 1);
prior = zeros(ncomb, 1);
comb = zeros(1, L);
for p = 0:ncomb-1
    tmp = p;
    for b = 1:L
        comb(b) = mod(tmp, nS) + 1;
        tmp = floor(tmp / nS);
    end
    if any(comb == 0)
        continue
    end
    all_comb{p+1} = comb;
    prior(p+1) = sum(comb ~= seq_i);  % distance from origin
end
all_comb = cat(1, all_comb{:});

all_seqs = cellstr(cell2mat(reshape(subs(all_comb), ncomb, [])));
all_seqs = strcat(flank{1}, all_seqs, flank{2});

[~, isort] = sort(prior, 'ascend');  % original seq first, followed by similar
all_seqs = all_seqs(isort);
