function resamp = resample_profiles(profiles, N)
% resamp = resample_profiles(profiles, N)
%   sample [N] reads-profiles based the given read distribution in
%   [profiles] (a cell array of gene profiles).
%   see also the corresponding function in RP.py.
%
% Alon Diament, Tuller Lab, July 2017.

reads_per_gene = cellfun(@nansum, profiles);
if abs(1 - sum(reads_per_gene) / N) < 1e-2
    fprintf('skipping resampling, sum reads = N\n');
    % avoid introducing noise if called for a dataset with its own sum of reads
    % (up to 1% error)
    resamp = profiles;
    return
end

% determine samples per gene (then, propagate to nucleotides)
% sampling is uniform, but non-deterministic
reads_per_gene = resample_reads(reads_per_gene, N);

nP = length(reads_per_gene);
resamp = cell(nP, 1);
err = 0;
for p = 1:nP
    nans = isnan(profiles{p});
    profiles{p}(nans) = 0;
    if sum(profiles{p}) > 0
        resamp{p} = resample_reads(profiles{p}, reads_per_gene(p));
    else
        resamp{p} = profiles{p};
    end
    resamp{p}(nans) = NaN;
    err = err + abs(nansum(resamp{p}) - reads_per_gene(p));
end
%    fprint('resample_profile: err = %.2f\n', err);

end
