function synCUFS = calc_synCUFS(seq)
% synCUFS = calc_synCUFS(seq)
%   computes the Synonymous Codon Usage Frequency Similarity (CUFS) for a given
%   cell array [seq] containing ORF sequences.
%   cite: 
%   Diament, Pinter, Tuller. Three-dimensional eukaryotic genomic 
%   organization is strongly correlated with codon usage expression and 
%   function. Nature Communications, 2014.
%
%   Alon Diament, Tuller Lab, 2013.

ng = length(seq);

codon_freq = zeros(ng, 64);
for i = 1:ng
    tmp = codonbias(seq{i});
    f = fieldnames(tmp);
    freq_vec = zeros(1, 64);
    ind(1) = 1;
    for j = 1:length(f)
        aa_freq = tmp.(f{j}).Freq;
        % replace NaNs (AA abscent) with uniform distribution
        aa_freq(isnan(aa_freq)) = 1/length(aa_freq);
        ind(2) = ind(1) + length(aa_freq) - 1;
        freq_vec(ind(1):ind(2)) = aa_freq;
        ind(1) = ind(2) + 1;
    end
    codon_freq(i, :) = freq_vec;
end

synCUFS = calc_ES_dist(codon_freq);
end


function D = calc_ES_dist(freqs)
% D = callc_ES_dist(freqs)
%   calculates the Endres-Schindelin all pair distances for the 
%   distributions in rows of matrix [freqs].

ng = size(freqs, 1);
D = zeros(ng, ng);

for i = 1:ng
    P = freqs(i, :);
    pi = P~=0;   % avoiding NaNs
    for j = 1:i-1
        Q = freqs(j, :);
        M = (P + Q)/2;
        qi = Q~=0;

        D(i, j) = sqrt( ...
                       sum( log( P(pi)./ M(pi)).* P(pi)) + ...
                       sum( log( Q(qi)./ M(qi)).* Q(qi)));
    end
end

D = D + D';
end