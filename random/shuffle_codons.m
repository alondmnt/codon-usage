function rseq_nt = shuffle_codons(seq_nt)
% rseq_nt = shuffle_codons(seq_nt)
%   shuffle codons within [seq_nt], preserving codon bias and AA seq.
%   when [seq_nt] is a cell array, shuffle between all genes.
%
% Alon Diament, Tuller Lab, 2016.

multi_gene = iscell(seq_nt);
if multi_gene
    lens = cellfun(@length, seq_nt);
    isbad = mod(lens, 3);
    if any(isbad)
        error('coding seq %s not divisable by 3', mat2str(find(isbad)));
    end
    seq_test = nt2aa(seq_nt, 'AlternativeStartCodons', false);
    seq_nt = [seq_nt{:}];
end

if mod(length(seq_nt), 3)
    error('coding seq not divisible by 3');
end

seq_aa = nt2aa(seq_nt, 'AlternativeStartCodons', false);
aa_list = fieldnames(aacount(seq_aa));

rseq_nt = seq_nt;
for a = 1:length(aa_list)
    pos_aa = find(seq_aa == aa_list{a});
    pos_nt = 1 + 3*(pos_aa-1);
    rpos_nt = pos_nt(randperm(length(pos_nt)));
    rseq_nt(pos_nt) = rseq_nt(rpos_nt);
    rseq_nt(pos_nt+1) = rseq_nt(rpos_nt+1);
    rseq_nt(pos_nt+2) = rseq_nt(rpos_nt+2);
end

if multi_gene
    lens = cumsum([1; lens(:)]);
    rseq_nt = arrayfun(@(i, j) rseq_nt(i:j-1), lens(1:end-1), lens(2:end), ...
                       'UniformOutput', false);
    if ~all(strcmp(nt2aa(rseq_nt), seq_test))
        error('too many cooks!');
    end
end
