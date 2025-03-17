function sources = strongest_sources(lf_brain, n_brain_sources, n)
    P2P = zeros(1, n_brain_sources);
    for i = 1:size(lf_brain, 2)
        P2P(:, i) = max(lf_brain(:, i)) - min(lf_brain(:, i));
    end
    P2P_sorted = sort(P2P, 'descend');
    strongest_sources = P2P_sorted(:, 1:n);
    sources = find(ismember(P2P, strongest_sources));
end

