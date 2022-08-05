function score_vec = calc_score_from_weights(w, seq)
% Alon Diament, Tuller Lab, 2015.

%% CALCULATION
if iscell(seq)
    % many genes
    score_vec = cellfun(@(x) one_score(w, x), seq);
else
    % single gene
    score_vec = one_score(w, seq);
end

end

function score = one_score(w, seq)

if iscell(seq)
    score = mean(cellfun(@(x) one_score(w, x), seq));
    return
end

gene_count = cell2mat(struct2cell(codoncount(seq)))';
nn = ~(isnan(w) | w==0);
score = exp(gene_count(nn) * log(w(nn)) / sum(gene_count(nn)));
% this computation of geometric mean was presented in Class 2 / Slide 22.

end
