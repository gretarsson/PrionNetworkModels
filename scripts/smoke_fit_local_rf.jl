#!/usr/bin/env julia

using PrionNetworkModels

function main()
    root = dirname(@__DIR__)
    n_samples = parse(Int, get(ENV, "PNM_SMOKE_SAMPLES", "150"))
    n_warmup = parse(Int, get(ENV, "PNM_SMOKE_WARMUP", "150"))
    spec = RunSpec(
        model = ModelSpec(name="LOCAL-RF", transport="retrograde", parameter_sharing="independent"),
        data = DataSpec(
            network = joinpath(root, "data/examples/network.csv"),
            observations = joinpath(root, "data/examples/observations.csv"),
            region_label_style = "synthetic",
        ),
        seeding = SeedingSpec(seed_indices=[1], infer_seed=false),
        inference = InferenceSpec(
            n_chains = 1,
            target_acceptance = 0.8,
            sampler = "NUTS",
            n_samples = n_samples,
            n_warmup = n_warmup,
            mean_data = false,
            ignore_seed = false,
        ),
        run_name = "smoke-fit-local-rf",
    )

    paths = fit_and_save_run(spec; run_root=joinpath(root, "runs"), run_id="smoke-fit-local-rf")
    println(paths.run_dir)
end

main()
