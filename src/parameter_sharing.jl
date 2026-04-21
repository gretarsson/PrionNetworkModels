function build_region_groups(labels::Vector{String})
    bases = map(l -> length(l) > 1 ? l[2:end] : l, labels)
    uniq = unique(bases)
    gid = Dict(base => i for (i, base) in enumerate(uniq))
    return [gid[b] for b in bases]
end
