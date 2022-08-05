function ENC = calc_ENC(seq)
% computes the Effective Number of Codons (Wright 1990).
%
% Alon Diament, Tuller Lab, August 2014.

if iscell(seq)
    ENC = cellfun(@calc_ENC, seq);
    return;
end

% using similar notations to (Wright, 1990)
P = codonbias(seq);
N = aacount(nt2aa(seq, 'ACGTOnly', false));
Fd = zeros(6,1); % vector for degeneracy types
Cd = zeros(6,1); % count degeneracy types

for f = fieldnames(P)'
    if strcmpi(f{1}, 'END')
        continue;
    end
    ThisP = P.(f{1});
    ThisN = N.(aminolookup(f{1}));
    if ThisN <= 1
        % div by zero
        continue;
    end
    ThisDeg = length(ThisP.Freq);
    % the mean over AA represents each degeneracy
    ThisF = (ThisN * sum(ThisP.Freq.^2) - 1) / (ThisN - 1);
    if ThisF == 0
        continue
    end
    Fd(ThisDeg) = Fd(ThisDeg) + ThisF;
    Cd(ThisDeg) = Cd(ThisDeg) + 1;
end

ind = Cd~=0;
Fd(ind) = Fd(ind) ./ Cd(ind);

% misssing AA cases
ind = find(Fd==0);
Fd(ind) = 1./ind;
if Cd(3) == 0
    Fd(3) = mean( Fd([2,4]));
end

ENC = 2 + 9/Fd(2) + 1/Fd(3) + 5/Fd(4) + 3/Fd(6);
ENC = min(61, ENC);
end
