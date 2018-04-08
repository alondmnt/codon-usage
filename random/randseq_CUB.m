function seqNT = randseq_CUB(seqAA, CUB)
% seqNT = rand_CUB(seqAA, CUB)
%  generate a random sequence with a given codon usage bias.
%  [CUB] can be a cellstr of gene sequences, or the output
%  of codonbias().
%
% Alon Diament, Tuller Lab, June 2017.

if ~isstruct(CUB)
    if all(~mod(cellfun(@length, CUB), 3))
        CUB = codonbias(strcat(CUB{:}));
    else
        error('too many cooks');
    end
end

nS = length(seqAA);
seqNT = cell(nS, 1);
for i = 1:nS
    seqNT{i} = rand_codon(CUB.(aminolookup(seqAA(i))));
end
seqNT = strcat(seqNT{:});
end


function codon = rand_codon(AA)
    cdf = cumsum(AA.Freq);
    codon = AA.Codon{find(rand <= cdf, 1)};
end
