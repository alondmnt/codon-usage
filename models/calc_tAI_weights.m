function [W, GCN] = calc_tAI_weights(tGCN, S, S_rules)
% [W, GCN] = calc_tAI_weights(tGCN, S, S_rules)
% model: tRNA Adaptation Index (dos Reis, 2004)
%
% tGCN: tRNA Gene Copy Numbers given according to anti-codons.
% S: s-values for all possible interactions between codons and tRNAs
% S_rules: all interaction rules, 1st anti-codon position (1st col) and 
%   3rd codon position (2nd col) and additional constraints on degree of AA
%   (3rd col).
% to compute the tAI score use: calc_score_from_weights()
%
% Alon Diament, Tuller Lab, 2015.

if nargin < 3
    S_rules = {'A',    'T',     1; ... % watson-crick
               'G',    'C',     1; ...
               'T',    'A',     1; ...
               'C',    'G',     1; ...
               'G',    'T',     2; ... % wobble
               'A',    'C',     3; ...
               'A',    'A',     3; ...
               'T',    'G',     2};
%                'C',    'A',     2};  % Lysidine
end
if nargin < 2
    S = [0, 0, 0, 0, 0.561, 0.28, 0.9999, 0.68, 0.89]; % 0.41 original / 0.561 optimization vs. PA
end

t_GCN = tGCN.GCN;
t_anti = tGCN.anti_codon;
t_codon = cellfun(@seqrcomplement, t_anti, ...
    'UniformOutput', false);

codon_list = fieldnames(codoncount(''));
stop_list = {'TGA', 'TAA', 'TAG'};
nC = length(codon_list);
W = zeros(nC, 1);

cod2aa = cellfun(@aa2int, nt2aa(codon_list, ...
                 'AlternativeStartCodons', false));
deg = arrayfun(@(x) nnz(cod2aa == x), cod2aa);

for c = 1:nC
    if any(strcmp(codon_list{c}, stop_list))
        W(c) = NaN;
        continue;
    end
    % candidate anti-codons (compare first 2 pos in codon)
    cand = find(strncmp(codon_list{c}, t_codon, 2))';
    if isempty(cand)
        continue;
    end
    % check S-rules for 3rd position
    for cand = cand
        is_tRNA_fits = @(x, y, z) (t_anti{cand}(1) == x) & ... % 1st anti pos
            (codon_list{c}(3) == y) & ... % 3rd codon pos
            (deg(c) >= z); % AA degree
        idx = cellfun(is_tRNA_fits, ... 
            S_rules(:, 1), S_rules(:, 2), S_rules(:, 3));
        if nnz(idx)
            W(c) = W(c) + (1 - S(idx))*t_GCN(cand);
        end
    end
end

W(strcmp('ATG', codon_list)) = t_GCN(strcmp('ATG', t_codon));
W = W / max(W);
W(W==0) = geomean(W(W~=0 & ~isnan(W)));

% additional output
GCN = zeros(length(W), 1);
[~, ind] = ismember(t_codon, codon_list);
GCN(ind) = t_GCN;

end
