function seq_nt = maximize_CUB(seq_nt, ref_nt)
% maximize_CUB(seq, ref)
%   replace all codons with optimal ones according to reference
%   sequences given in [ren_nt] (cell array). alternatively,
%   [ref_nt] can also be the struct output of codonbias().
%
% Alon Diament, Tuller Lab, 2016.

seq_aa = nt2aa(seq_nt, 'AlternativeStartCodons', false);

if isstruct(ref_nt)
    CUB = ref_nt;
else
    lens = cellfun(@length, ref_nt);
    ind = ~mod(lens, 3);
    fprintf('maximize_CUB: ignored %d ref seqs\n', sum(~ind));
    CUB = codonbias(strcat(ref_nt{ind}));  % ensuring division by 3
end

aa_list = fieldnames(CUB);
for a = 1:length(aa_list)
    aa = aminolookup(aa_list{a});
    pos = 3*(find(seq_aa == aa) - 1) + 1;  % aa to codon
    [~, bestcode] = max(CUB.(aa_list{a}).Freq);
    for p = pos
        seq_nt(p:p+2) = CUB.(aa_list{a}).Codon{bestcode};
    end
end

% test
if ~all(nt2aa(seq_nt, 'AlternativeStartCodons', false) == seq_aa)
    error('too many cooks!');
end
