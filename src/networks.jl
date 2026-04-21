function laplacian_out(W::AbstractMatrix; self_loops::Bool=false, retro::Bool=false)
    Wwork = Matrix{Float64}(W)
    if retro
        Wwork = transpose(Wwork)
    end
    N = size(Wwork, 1)
    if !self_loops
        for i in 1:N
            Wwork[i, i] = 0.0
        end
    end
    D = Diagonal(vec(sum(Wwork, dims=2)))
    return D - Wwork
end

function read_network_csv(path::AbstractString)
    raw = readdlm(path, ',')
    labels = String.(raw[1, 2:end])
    W = Matrix{Float64}(raw[2:end, 2:end])
    return (W = W, labels = labels)
end

function build_transport_operator(path::AbstractString; transport::AbstractString="retrograde", self_loops::Bool=false)
    network = read_network_csv(path)
    W = copy(network.W)
    if any(W .> 0)
        W ./= maximum(W[W .> 0])
    end
    mode = String(transport)
    if mode == "retrograde"
        L = transpose(laplacian_out(W; self_loops=self_loops, retro=true))
    elseif mode == "anterograde"
        L = transpose(laplacian_out(W; self_loops=self_loops, retro=false))
    elseif mode == "euclidean"
        L = transpose(laplacian_out(W; self_loops=self_loops, retro=true))
    elseif mode == "bidirectional"
        Lr = transpose(laplacian_out(W; self_loops=self_loops, retro=true))
        La = transpose(laplacian_out(W; self_loops=self_loops, retro=false))
        L = (Lr + La) ./ 2
    else
        error("Unsupported transport mode: $transport")
    end
    return (L = Matrix{Float64}(L), labels = network.labels, W = W)
end
