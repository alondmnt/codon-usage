function resamp = resample_reads(read_distrib, N)
% resamp = resample_reads(read_distrib, N)
%   resample [N] reads from the given distribution [read_distrib].
%   used in resample_profiles()
%   see also the corresponding function in RP.py.
%
% Alon Diament, Tuller Lab, July 2017.

if size(read_distrib, 1) < size(read_distrib, 2)
    read_distrib = read_distrib';
end

M = length(read_distrib);
% quick inverse transform sampling using sorting
cdf = cumsum(read_distrib) / sum(read_distrib);
[~, I] = sort([0; cdf; rand(N, 1)]);
resamp(I, 1) = 1:(M + N + 1);  % get the ranks
resamp = resamp(1:M+1);  % of cdf values
% cool trick for counting the no. of drawn reads between cdf values
resamp = diff(resamp) - 1;  % (index fix, o/w sum != N)

end
