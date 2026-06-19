#!/usr/bin/env julia

using CSV
using DataFrames
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RUN_ROOT = joinpath(PROJECT_ROOT, "runs")
const OUT_DIR = joinpath(PROJECT_ROOT, "paper-rf", "results", "model_figures")

function load_run_waic(run_id::AbstractString)
    run_dir = joinpath(RUN_ROOT, run_id)
    summary_path = joinpath(run_dir, "waic_summary.csv")
    pointwise_path = joinpath(run_dir, "waic_pointwise.csv")
    isfile(summary_path) || error("Missing WAIC summary for run: $run_id")
    isfile(pointwise_path) || error("Missing pointwise WAIC table for run: $run_id")
    summary = CSV.read(summary_path, DataFrame)
    pointwise = CSV.read(pointwise_path, DataFrame)
    return (
        waic = Float64(summary.waic[1]),
        waic_i = Float64.(pointwise.waic_i),
    )
end

function paired_delta_waic(rowspecs)
    loaded = [load_run_waic(spec.run_id) for spec in rowspecs]
    waic = [item.waic for item in loaded]
    waic_i = [item.waic_i for item in loaded]
    best_idx = argmin(waic)
    ref = waic_i[best_idx]
    rows = DataFrame()
    for (idx, spec) in enumerate(rowspecs)
        d_i = waic_i[idx] .- ref
        mask = .!(isnan.(d_i) .| isinf.(d_i))
        delta = idx == best_idx ? 0.0 : sum(d_i[mask])
        se = idx == best_idx ? 0.0 : sqrt(length(d_i[mask]) * var(d_i[mask]))
        class = idx == best_idx ? "best" : ((delta - 2se) <= 0 <= (delta + 2se) ? "tied" : "worse")
        push!(rows, (
            panel = spec.panel,
            model = spec.model,
            run_id = spec.run_id,
            waic = waic[idx],
            delta_waic = delta,
            se_delta_waic = se,
            class = class,
        ))
    end
    return rows
end

function export_transport_waic()
    specs = [
        [
            (panel = "DIFF", model = "euclidean", run_id = "striatum_DIFF_EUCL"),
            (panel = "DIFF", model = "anterograde", run_id = "striatum_DIFF_ANTERO"),
            (panel = "DIFF", model = "retrograde", run_id = "striatum_DIFF_RETRO_paper"),
            (panel = "DIFF", model = "bidirectional", run_id = "striatum_DIFF_BIDIR"),
        ],
        [
            (panel = "DIFF-R", model = "euclidean", run_id = "striatum_DIFF-R_EUCL"),
            (panel = "DIFF-R", model = "anterograde", run_id = "striatum_DIFF-R_ANTERO"),
            (panel = "DIFF-R", model = "retrograde", run_id = "striatum_DIFF-R_RETRO_paper"),
            (panel = "DIFF-R", model = "bidirectional", run_id = "striatum_DIFF-R_BIDIR"),
        ],
        [
            (panel = "DIFF-RF", model = "euclidean", run_id = "striatum_DIFF-RF_EUCL"),
            (panel = "DIFF-RF", model = "anterograde", run_id = "striatum_DIFF-RF_ANTERO"),
            (panel = "DIFF-RF", model = "retrograde", run_id = "striatum_DIFF-RF_RETRO_paper"),
            (panel = "DIFF-RF", model = "bidirectional", run_id = "striatum_DIFF-RF_BIDIR"),
        ],
    ]
    CSV.write(joinpath(OUT_DIR, "figure2_transport_waic.csv"), reduce(vcat, paired_delta_waic.(specs)))
end

function export_model_class_waic()
    specs = [
        (panel = "striatum", model = "DIFF", run_id = "striatum_DIFF_RETRO_paper"),
        (panel = "striatum", model = "DIFF-R", run_id = "striatum_DIFF-R_RETRO_paper"),
        (panel = "striatum", model = "DIFF-RF", run_id = "striatum_DIFF-RF_RETRO_paper"),
    ]
    CSV.write(joinpath(OUT_DIR, "figure3_model_waic.csv"), paired_delta_waic(specs))
end

function export_null_waic()
    connectivity = CSV.read(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_connectivity_nulls", "waic_values.csv"), DataFrame)
    seeding = CSV.read(joinpath(RUN_ROOT, "striatum_DIFF-RF_RETRO_seed_nulls", "waic_values.csv"), DataFrame)
    empirical = load_run_waic("striatum_DIFF-RF_RETRO_paper").waic
    rows = DataFrame(null_type=String[], source=String[], waic=Float64[], is_empirical=Bool[])
    for row in eachrow(connectivity)
        push!(rows, ("connectivity", String(row.source), Float64(row.waic), false))
    end
    for row in eachrow(seeding)
        push!(rows, ("seeding", String(row.source), Float64(row.waic), false))
    end
    push!(rows, ("connectivity", "striatum_DIFF-RF_RETRO_paper", empirical, true))
    push!(rows, ("seeding", "striatum_DIFF-RF_RETRO_paper", empirical, true))
    CSV.write(joinpath(OUT_DIR, "figure4_null_waic.csv"), rows)
end

function main()
    mkpath(OUT_DIR)
    export_transport_waic()
    export_model_class_waic()
    export_null_waic()
    println("Wrote WAIC tables from run bundles to $OUT_DIR")
end

main()
