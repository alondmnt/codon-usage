function RCBS = calc_RCBS(seq)
% compute the Relative Codon Bias Score (Roymondal, 2009)
%
% Alon Diament, Tuller Lab, August 2014.

if iscell(seq)
    RCBS = cellfun(@calc_RCBS, seq);
    return;
end

CodonCount = cell2mat(struct2cell(codoncount(seq))); % sorted
P = CodonCount / sum(CodonCount);
BCC = calc_BCC(calc_BNC(seq)); % background codon composition
D = (P - BCC) ./ BCC; % per-codon score

nz = P~=0;
RCBS = log(1 + D(nz)') * CodonCount(nz) / sum(CodonCount);
RCBS = exp(RCBS) - 1; % gene score
end


function [BNC, NTs] = calc_BNC(back_seq)
% background nucleotide composition (BNC)

BNC = zeros(4,3); % 4-NTs x 3-positions
NTs = {'A','C','G','T'};

for pos = 1:3
    seq = back_seq(pos:3:end);   
    % derive NT distibution directly
    for nt = 1:4
        BNC(nt, pos) = sum(seq == NTs{nt}) / length(seq);
    end
end
end


function BCC = calc_BCC(BNC)
% background codon composition (BCC)
% given the background nucleotide composition (BNC)

BCC = zeros(64,1);
i = 0;

for n1 = 1:4
    f(1) = BNC(n1,1);
    for n2 = 1:4
        f(2) = BNC(n2,2);
        for n3 = 1:4
            f(3) = BNC(n3,3);
            % codon probability is f(1)*f(2)*f(3)
            i = i + 1;
            BCC(i) = f(1)*f(2)*f(3);
        end
    end
end
end
