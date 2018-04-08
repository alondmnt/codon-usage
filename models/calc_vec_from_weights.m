function vec = calc_vec_from_weights(seq, w, codon_list)
% vec = calc_vec_from_weights(seq, w, codon_list)
%  generate an array of values for a given gene according to the
%  given codon-specific weights.
%
% Alon Diament, Tuller Lab, May 2017.

assert(all(ismember(seq, 'AGCTacgt')), 'bad sequence');

if nargin < 3
    codon_list = fieldnames(codoncount(''));  % lexicographic order
end

codons = cellstr(reshape(seq, 3, [])');
[~, ind] = ismember(upper(codons), codon_list);
vec = w(ind(ind > 0))';

end
