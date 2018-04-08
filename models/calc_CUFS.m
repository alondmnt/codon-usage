function CUFS = calc_CUFS(seq)
% CUFS = calc_CUFS(seq)
%   computes the Codon Usage Frequency Similarity (CUFS) for a given
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
    [~, tmp] = codoncount(seq{i});
    codon_freq(i, :) = reshape(tmp, 1, 64) / sum(tmp(:));
end

CUFS = calc_ES_dist(codon_freq);
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
