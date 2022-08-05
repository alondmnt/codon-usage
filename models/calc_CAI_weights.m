function w = calc_CAI_weights(ref)
% w = calc_CAI_weights(ref)
% model: Codon Adaptation Index (Sharp & Li, 1987)
%
%   computes the CAI weights for based on the reference set of sequences [ref].
% returns:
%   w:      CAI weights vector for the given reference set.
%
% to compute the CAI score use: calc_score_from_weights()
%
% Alon Diament, Tuller Lab, 2015.

%% PROCESS REFERENCE
ref_seq = [];
for i = 1:length(ref)
    last_legal_codon = length(ref{i});
    last_legal_codon = last_legal_codon - mod(last_legal_codon, 3); % length is now divisible by 3
    ref_seq = strcat(ref_seq, ref{i}(1:last_legal_codon));
end

ref_count = codoncount(ref_seq);
codon_list = fieldnames(ref_count);
cod2aa = cellfun(@aa2int, nt2aa(codon_list, ...
    'AlternativeStartCodons', false));
ref_count = cell2mat(struct2cell(ref_count))';
ref_count(ref_count == 0) = 0.5; % in case a codon is missing
 
%% NORMALIZED WEIGHTS
w = nan(length(codon_list), 1);
for i = 1:20
    % find all codons of AA(i)
    idx = cod2aa == i;
    w(idx) = ref_count(idx) / max(ref_count(idx));
end
% note that the array [w] has the same order of values as [ref_count] and
% later [gene_count]. also note, that stop codons were ignored and we have
% NaN's left in their positions.
 
end
