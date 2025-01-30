"""
"""

"""
    IcaWorkspace

Structure containing pre-allocated arrays for ICA computation.
"""
struct IcaWorkspace
    # Basic arrays
    weights::Matrix{Float64}
    oldweights::Matrix{Float64}
    startweights::Matrix{Float64}
    bias::Vector{Float64}
    data_block::Matrix{Float64}
    u::Matrix{Float64}
    y::Matrix{Float64}
    delta_weights::Matrix{Float64}
    temp_matrix::Matrix{Float64}
    yu::Matrix{Float64}
    uu::Matrix{Float64}
    BI::Matrix{Float64}
    onesrow::Matrix{Float64}
    permute_indices::Vector{Int}
    
    # Weight change tracking
    olddelta::Vector{Float64}
    degconst::Float64
    
    # Extended ICA specific
    signs::Union{Nothing, Diagonal{Float64}}
    old_kurt::Union{Nothing, Vector{Float64}}
    oldsigns::Union{Nothing, Diagonal{Float64}}
    partact::Union{Nothing, Matrix{Float64}}
    kk::Union{Nothing, Vector{Float64}}
    m2::Union{Nothing, Vector{Float64}}
    m4::Union{Nothing, Vector{Float64}}
    new_signs::Union{Nothing, Vector{Float64}}

    function IcaWorkspace(n_channels::Int, n_samples::Int, block::Int, extended::Bool, ext_params::Union{Nothing, ExtendedIcaPrms})
        # Basic arrays
        ws = new(
            Matrix{Float64}(I, n_channels, n_channels),  # weights
            similar(Matrix{Float64}(I, n_channels, n_channels)),  # oldweights
            similar(Matrix{Float64}(I, n_channels, n_channels)),  # startweights
            zeros(n_channels),  # bias
            zeros(n_channels, block),  # data_block
            zeros(n_channels, block),  # u
            zeros(n_channels, block),  # y
            zeros(n_channels, n_channels),  # delta_weights
            zeros(n_channels, n_channels),  # temp_matrix
            zeros(n_channels, n_channels),  # yu
            zeros(n_channels, n_channels),  # uu
            block * Matrix{Float64}(I, n_channels, n_channels),  # BI
            ones(1, block),  # onesrow
            Vector{Int}(undef, n_samples),  # permute_indices
            zeros(n_channels * n_channels),  # olddelta
            180.0 / π,  # degconst
            nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing  # Extended ICA fields
        )

        # Copy initial weights
        copyto!(ws.oldweights, ws.weights)
        copyto!(ws.startweights, ws.weights)

        # Initialize extended ICA arrays if needed
        if extended
            ws.signs = Diagonal(vcat(-ones(ext_params.n_subgauss), ones(n_channels - ext_params.n_subgauss)))
            ws.old_kurt = zeros(n_channels)
            ws.oldsigns = copy(ws.signs)
            ws.partact = zeros(n_channels, min(ext_params.kurt_size, n_samples))
            ws.kk = zeros(n_channels)
            ws.m2 = zeros(n_channels)
            ws.m4 = zeros(n_channels)
            ws.new_signs = zeros(n_channels)
        end

        return ws
    end
end

function infomax_ica(data::Matrix{Float64};
    extended::Bool = false,
    n_components::Union{Nothing,Int} = nothing,
    max_iter::Int = 512,
    l_rate::Union{Nothing,Float64} = nothing,
    w_change::Float64 = 1e-12,
    anneal_deg::Float64 = 60.0,
    anneal_step::Float64 = 0.9,
    use_bias::Bool = true,
    do_whitening::Bool = true,
    verbose::Bool = true)

    # Input validation
    n_channels, n_samples = size(data)
    if !isnothing(n_components)
        if n_components > n_channels
            throw(ArgumentError("n_components ($n_components) cannot exceed number of channels ($n_channels)"))
        end
        if n_components < 1
            throw(ArgumentError("n_components must be positive"))
        end
    else
        n_components = n_channels
    end

    # Set default learning rate if not provided
    if isnothing(l_rate)
        l_rate = min(0.00065/log(n_channels), 0.1)
    end

    # # Configure logging
    # if verbose
    #     logger = ConsoleLogger(stdout, Logging.Info)
    # else
    #     logger = ConsoleLogger(stdout, Logging.Warn)
    # end
    # global_logger(logger)

    # Initialize parameters
    params = IcaPrms(
        l_rate = l_rate,
        max_iter = max_iter,
        w_change = w_change,
        use_bias = use_bias,
        anneal_deg = anneal_deg,
        anneal_step = anneal_step
    )
    
    ext_params = extended ? ExtendedIcaPrms() : nothing

    # Preprocessing
    @info "Preprocessing data..."
    data_copy = copy(data)  # Preserve original data
    
    # PCA reduction if requested
    if n_components < n_channels
        @info "Reducing dimensions from $n_channels to $n_components using PCA"
        data_reduced, pca_matrix = pca_reduction(data_copy, n_components)
    else
        data_reduced = data_copy
        pca_matrix = nothing
    end

    # Whitening
    if do_whitening
        @info "Whitening data..."
        data_whitened, sphere = pre_whiten(data_reduced, return_scale=true)
    else
        data_whitened = data_reduced
        sphere = nothing
    end

    # Initialize workspace
    @info "Initializing ICA..."
    workspace = initialize_workspace(data_whitened, extended, ext_params)

    # Run ICA iterations
    @info "Running ICA iterations..."
    weights = run_ica_iterations!(data_whitened, workspace, params, ext_params)

    # Post-processing
    @info "Post-processing results..."
    
    # Calculate mixing matrix (topography)
    topo = pinv(weights)

    # If PCA was applied, adjust matrices
    if !isnothing(pca_matrix)
        topo = pca_matrix' * topo
        weights = weights * pca_matrix
    end

    # If whitening was applied, adjust matrices
    if !isnothing(sphere)
        topo = topo * sphere
        weights = inv(sphere) * weights
    end

    # Generate labels
    n_out = size(weights, 1)
    labels = ["IC$i" for i = 1:n_out]

    @info "ICA completed successfully"
    
    # Return results with additional matrices
    return InfoIca(
        topo,           # Mixing matrix
        weights,        # Unmixing matrix
        labels,         # Component labels
        sphere,         # Sphering matrix (if used)
        pca_matrix      # PCA matrix (if used)
    )
end