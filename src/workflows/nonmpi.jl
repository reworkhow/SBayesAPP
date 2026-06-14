using DelimitedFiles
using Distributions
using LinearAlgebra
using Base.Threads
using ProgressMeter
using Random
using Statistics
using Dates
using JLD2

function make_nonmpi_block_workspace(blocks, block_range, nCategory::Int, nlabel::Int, nTraits::Int)
    max_nMarker = maximum(size(blocks[b].x_arrays[end], 2) for b in block_range)
    max_nEigen = maximum(size(blocks[b].x_arrays[end], 1) for b in block_range)
    return (
        local_nLoci_array_vec=[fill(0, nlabel) for _ in 1:nCategory],
        beta=zeros(nTraits),
        new_alpha=zeros(nTraits),
        old_alpha=zeros(nTraits),
        residual=zeros(nTraits),
        delta=zeros(nTraits),
        marker_order=collect(1:max_nMarker),
        what_cat_1=zeros(max_nEigen),
        what_cat_2=zeros(max_nEigen),
        what_total_1=zeros(max_nEigen),
        what_total_2=zeros(max_nEigen),
        residual_sse=zeros(nTraits, nTraits),
    )
end

function reset_nonmpi_block_workspace!(workspace)
    for counts in workspace.local_nLoci_array_vec
        fill!(counts, 0)
    end
    return workspace
end

function accumulate_matrices!(dest, matrices)
    fill!(dest, 0.0)
    for matrix in matrices
        @inbounds for index in eachindex(dest, matrix)
            dest[index] += matrix[index]
        end
    end
    return dest
end

function accumulate_vectors!(dest, vectors)
    fill!(dest, 0.0)
    for vector in vectors
        @inbounds for index in eachindex(dest, vector)
            dest[index] += vector[index]
        end
    end
    return dest
end

function run_nonmpi_block_range!(
    blocks,
    block_range,
    betaArray,
    alphaArray,
    deltaArray,
    meanAlpha,
    Ainv_vec,
    Pi,
    marker_probit_tree_state,
    R_blk,
    varg_blk_cat,
    varg_cov_blk_cat,
    totalvarg_blk,
    ssq_blk_cat,
    workspace;
    annotation_prior_model,
    effect_nCon,
    my_nsnp,
    nCategory,
    nTraits,
    nlabel,
    do_thin,
    iIter,
    need_block_totalvarg,
    estimate_vare,
    df_R,
    scale_R,
)
    nTraits == 2 || error("run_nonmpi_block_range! currently expects exactly 2 traits.")

    reset_nonmpi_block_workspace!(workspace)
    local_nLoci_array_vec = workspace.local_nLoci_array_vec
    β = workspace.beta
    newα = workspace.new_alpha
    oldα = workspace.old_alpha
    w = workspace.residual
    δ = workspace.delta
    marker_order = workspace.marker_order
    what_cat_1 = workspace.what_cat_1
    what_cat_2 = workspace.what_cat_2
    what_total_1 = workspace.what_total_1
    what_total_2 = workspace.what_total_2
    residual_sse = workspace.residual_sse

    for b in block_range
        block = blocks[b]
        xpx_vec = block.xpx::Vector{Vector{Float64}}
        xArray_vec = block.x_arrays::Vector{Matrix{Float64}}
        wArray = block.transformed_y
        category_indices_by_marker = block.annotation_category_indices::Vector{Vector{Int}}
        category_values_by_marker = block.annotation_category_values::Vector{Vector{Float64}}
        category_effect_indices_by_marker = block.annotation_effect_indices::Vector{Vector{Int}}
        category_marker_indices_by_cat = block.category_marker_indices::Vector{Vector{Int}}
        category_effect_indices_by_cat = block.category_effect_indices::Vector{Vector{Int}}
        SNPIndexb = block.snp_indices
        nEigenb, nMarkerb = size(xArray_vec[end])
        nInd = block.n_gwas

        # initialize xArrayc/xpxc and block-specific R matrix 
        xArrayc = xArray_vec[1]
        xpxc = xpx_vec[1]

        r11 = R_blk[b][1, 1]
        r12 = R_blk[b][1, 2]
        r21 = R_blk[b][2, 1]
        r22 = R_blk[b][2, 2]
        rdet = r11 * r22 - r12 * r21 # determinant of R_blk[b]
        Rinv11 = r22 / rdet
        Rinv12 = -r12 / rdet
        Rinv21 = -r21 / rdet
        Rinv22 = r11 / rdet

        resize!(marker_order, nMarkerb)
        for marker in 1:nMarkerb
            marker_order[marker] = marker
        end
        shuffle!(marker_order)

        for marker in marker_order
            true_marker_num = SNPIndexb[marker]
            marker_category_indices = category_indices_by_marker[marker]
            marker_category_values = category_values_by_marker[marker]
            marker_effect_indices = category_effect_indices_by_marker[marker]
            for category_position in eachindex(marker_category_indices)
                cat = marker_category_indices[category_position]
                annotation_value = marker_category_values[category_position]
                annotation_value == 0.0 && continue
                design_index = (effect_nCon != 0 && cat <= effect_nCon + 1) ? cat : length(xArray_vec)
                xArrayc = xArray_vec[design_index]
                xpxc = xpx_vec[design_index]
                Ginv = Ainv_vec[cat]
                effect_index = marker_effect_indices[category_position]
                x = view(xArrayc, :, marker)
                xpx_marker = xpxc[marker]
                for trait = 1:nTraits
                    β[trait] = betaArray[trait][cat][effect_index]
                    oldα[trait] = newα[trait] = alphaArray[trait][cat][effect_index]
                    δ[trait] = deltaArray[trait][cat][effect_index]
                    w[trait] = dot(x, wArray[trait]) + xpx_marker * oldα[trait]
                end

                for k = 1:nTraits
                    other = k == 1 ? 2 : 1
                    Ginv11 = Ginv[k, k]
                    Ginv12 = Ginv[k, other]
                    Rinvkk = k == 1 ? Rinv11 : Rinv22
                    Rinvkother = k == 1 ? Rinv12 : Rinv21
                    C11 = Ginv11 + Rinvkk * xpx_marker
                    C12 = Ginv12 + xpx_marker * δ[other] * Rinvkother

                    # compute the conditional distribution parameters for the two possible states of δ[k]
                    # when δ[k] = 0, the effect is 0 and does not contribute to the likelihood, so the conditional distribution is determined by the prior and the contribution of the other trait (if it is nonzero)
                    invLhs0 = 1 / Ginv11
                    rhs0 = -Ginv12 * β[other]
                    gHat0 = rhs0 * invLhs0
                    # when δ[k] = 1, the effect is nonzero and contributes to the likelihood, so the conditional distribution is determined by both the prior and the likelihood contribution of this trait and the other trait (if it is nonzero)
                    invLhs1 = 1 / C11
                    rhs1 = (k == 1 ? (w[1] * Rinv11 + w[2] * Rinv21) : (w[1] * Rinv12 + w[2] * Rinv22)) - C12 * β[other]
                    gHat1 = rhs1 * invLhs1

                    d0 = k == 1 ? (0.0, δ[2]) : (δ[1], 0.0)
                    d1 = k == 1 ? (1.0, δ[2]) : (δ[1], 1.0)

                    logPrior0 = log_marker_state_prior(annotation_prior_model, Pi, marker_probit_tree_state, cat, true_marker_num, d0)
                    logPrior1 = log_marker_state_prior(annotation_prior_model, Pi, marker_probit_tree_state, cat, true_marker_num, d1)
                    logDelta0 = -0.5 * (log(Ginv11) - gHat0^2 * Ginv11) + logPrior0
                    logDelta1 = -0.5 * (log(C11) - gHat1^2 * C11) + logPrior1
                    probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1))

                    # sample δ[k] from its conditional distribution, then sample the effect size if δ[k] = 1, and update w accordingly

                    if rand() < probDelta1
                        δ[k] = 1
                        β[k] = newα[k] = gHat1 + randn() * sqrt(invLhs1)
                        LinearAlgebra.axpy!(oldα[k] - newα[k], x, wArray[k])
                    else
                        β[k] = gHat0 + randn() * sqrt(invLhs0)
                        δ[k] = 0
                        newα[k] = 0
                        if oldα[k] != 0
                            LinearAlgebra.axpy!(oldα[k], x, wArray[k])
                        end
                    end
                end

                local_nLoci_array_vec[cat][two_trait_state_index(δ[1], δ[2])] += 1

                for trait = 1:nTraits
                    betaArray[trait][cat][effect_index] = β[trait]
                    deltaArray[trait][cat][effect_index] = δ[trait]
                    alphaArray[trait][cat][effect_index] = newα[trait]
                    if do_thin
                        meanAlpha[trait][cat][effect_index] += (newα[trait] - meanAlpha[trait][cat][effect_index]) * iIter
                    end
                end
            end # for category_position in eachindex(marker_category_indices)
        end # for marker in marker_order
        block.transformed_y = wArray

        if need_block_totalvarg
            what_total_1_block = view(what_total_1, 1:nEigenb)
            what_total_2_block = view(what_total_2, 1:nEigenb)
            fill!(what_total_1_block, 0.0)
            fill!(what_total_2_block, 0.0)
            what_cat_1_block = view(what_cat_1, 1:nEigenb)
            what_cat_2_block = view(what_cat_2, 1:nEigenb)

            for cat in 1:nCategory
                XAb = (effect_nCon != 0 && cat <= effect_nCon + 1) ? xArray_vec[cat] : xArray_vec[end]
                fill!(what_cat_1_block, 0.0)
                fill!(what_cat_2_block, 0.0)
                cat_marker_indices = category_marker_indices_by_cat[cat]
                cat_effect_indices = category_effect_indices_by_cat[cat]
                ssq_trait_1 = 0.0
                ssq_trait_2 = 0.0
                for position in eachindex(cat_marker_indices)
                    marker = cat_marker_indices[position]
                    effect_index = cat_effect_indices[position]
                    alpha_trait_1 = alphaArray[1][cat][effect_index]
                    alpha_trait_2 = alphaArray[2][cat][effect_index]
                    ssq_trait_1 += alpha_trait_1^2
                    ssq_trait_2 += alpha_trait_2^2
                    alpha_trait_1 != 0.0 && LinearAlgebra.axpy!(alpha_trait_1, view(XAb, :, marker), what_cat_1_block)
                    alpha_trait_2 != 0.0 && LinearAlgebra.axpy!(alpha_trait_2, view(XAb, :, marker), what_cat_2_block)
                end

                if do_thin
                    varg_blk_cat[b][1, cat] = dot(what_cat_1_block, what_cat_1_block)
                    varg_blk_cat[b][2, cat] = dot(what_cat_2_block, what_cat_2_block)
                    varg_cov_blk_cat[b][cat] = dot(what_cat_1_block, what_cat_2_block)
                end
                if estimate_vare
                    ssq_blk_cat[b][1, cat] = ssq_trait_1
                    ssq_blk_cat[b][2, cat] = ssq_trait_2
                end
                LinearAlgebra.axpy!(1.0, what_cat_1_block, what_total_1_block)
                LinearAlgebra.axpy!(1.0, what_cat_2_block, what_total_2_block)
            end

            totalvarg_blk[b][1, 1] = dot(what_total_1_block, what_total_1_block)
            totalvarg_blk[b][2, 2] = dot(what_total_2_block, what_total_2_block)
            totalvarg_blk[b][1, 2] = totalvarg_blk[b][2, 1] = dot(what_total_1_block, what_total_2_block)
        end

        if estimate_vare
            sampled_R = sample_variance_sumstats!(residual_sse, wArray, nEigenb, df_R, scale_R)
            R_blk[b] = sampled_R
            Rcor = compute_correlation(sampled_R)
            for traiti = 1:nTraits
                thres = sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]
                if thres > 1.1
                    R_blk[b][traiti, traiti] = sampled_R[traiti, traiti]
                else
                    R_blk[b][traiti, traiti] = 1.0 / nInd[traiti]
                end
            end
            Rcov = Rcor * sqrt(R_blk[b][1, 1] * R_blk[b][2, 2])
            R_blk[b][1, 2] = R_blk[b][2, 1] = Rcov
        end
    end # for b in block_range

    return local_nLoci_array_vec
end

function build_nonmpi_run_context(config::ConfigTypes.NonMPIConfig)
    settings = build_nonmpi_sampler_settings(config)
    annotation_metadata = load_annotation_metadata(config.data_path, config.annot_file; nCon=settings.effective_n_con)
    start_pi_result = build_start_pi(config.st_path; estimate_pi=settings.estimate_pi)
    gprior_vec = config.annotation_prior_model == :group_dirichlet ?
        build_gprior_vec(config.st_path, annotation_metadata.nLoci_annot, start_pi_result.Pi00) :
        nothing

    return (
        config=config,
        annotation_metadata=annotation_metadata,
        settings=settings,
        startPi=start_pi_result.startPi,
        startPi00=start_pi_result.Pi00,
        Gprior_vec=gprior_vec,
    )
end

function write_nonmpi_run_config(analysis_path::AbstractString, context)
    config = context.config
    settings = context.settings
    run_config = [
        "created_at" => string(now()),
        "data_path" => config.data_path,
        "analysis_path" => config.analysis_path,
        "annot_file" => config.annot_file,
        "annot_dict" => config.annot_dict,
        "starting_value_dir" => config.starting_value_dir,
        "gscale_value_dir" => config.gscale_value_dir,
        "st_path" => config.st_path,
        "n_iter" => config.nIter,
        "burnin" => config.burnin,
        "thin" => config.thin,
        "seed" => config.seed,
        "n1" => config.n1,
        "n2" => config.n2,
        "effective_n_con" => settings.effective_n_con,
        "annotation_prior_model" => String(config.annotation_prior_model),
        "is_continue" => settings.is_continue,
        "estimate_vare" => settings.estimate_vare,
        "estimate_vara" => settings.estimate_vara,
        "estimate_pi" => settings.estimate_pi,
        "effective_estimate_gscale" => settings.estimate_Gscale,
        "estgscale_iter" => settings.estGscale_iter,
        "report_pleiotropic_qtl_effect_matrix" => config.report_pleiotropic_qtl_effect_matrix,
        "output_mcmc_delta" => config.output_mcmc_delta,
        "julia_threads" => nthreads(),
    ]
    open(joinpath(analysis_path, "run_config.toml"), "w") do io
        for (key, value) in run_config
            println(io, key, " = ", value)
        end
    end
end

function run_nonmpi_sampler!(context)
    config = context.config
    annotation_metadata = context.annotation_metadata
    settings = context.settings
    Gprior_vec = context.Gprior_vec

    analysis_path = config.analysis_path
    nIter = config.nIter
    starting_value_dir = config.starting_value_dir
    thin = config.thin
    is_continue = settings.is_continue
    annotation_prior_model = config.annotation_prior_model

    nCon = annotation_metadata.nCon
    estimate_vare = settings.estimate_vare
    estimate_vara = settings.estimate_vara
    estimate_pi = settings.estimate_pi
    estimate_Gscale = settings.estimate_Gscale
    estGscale_iter = settings.estGscale_iter
    report_pleiotropic_qtl_effect_matrix = config.report_pleiotropic_qtl_effect_matrix
    output_mcmc_delta = config.output_mcmc_delta

    ############################################################################
    #read data in current rank
    ############################################################################
    my_rank = 0

    block_data = load_nonmpi_block_data(config.data_path, config.annot_dict)
    blocks = block_data.blocks
    my_nblk = block_data.nblk
    my_nsnp = block_data.nsnp

    burnin = config.burnin

    nCategory = annotation_prior_model == :group_dirichlet ? annotation_metadata.nCat + nCon : 1
    nLoci_annot = annotation_prior_model == :group_dirichlet ? annotation_metadata.nLoci_annot : Int[my_nsnp]
    effect_nCon = annotation_prior_model == :group_dirichlet ? nCon : 0
    if Gprior_vec === nothing
        Gprior_vec = build_gprior_vec(config.st_path, nLoci_annot, context.startPi00)
    end
    
    nTraits = 2
    parameter_state = initialize_nonmpi_parameter_state(
        my_rank,
        nCategory,
        nTraits,
        Gprior_vec,
        context.startPi;
        n1=config.n1,
        n2=config.n2,
        estimate_vare=estimate_vare,
        estimate_vara=estimate_vara,
        estimate_Gscale=estimate_Gscale,
        is_continue=is_continue,
        starting_value_dir=starting_value_dir,
        gscale_value_dir=config.gscale_value_dir,
    )
    (; A_vec, Ainv_vec, Pi, Rprior, df_R, scale_R, df_G, scale_G_vec, estimate_Gscale) = parameter_state

    mkpath(analysis_path)

    
    write_nonmpi_run_config(analysis_path, context)
    writedlm(analysis_path * "annotationName.txt", annotation_metadata.annotationName)

    marker_probit_tree_state = nothing
    if annotation_prior_model == :group_dirichlet
        prepare_block_state!(blocks, effect_nCon, annotation_metadata.annotationType)
    else
        annotation_design = prepare_marker_probit_tree_block_state!(blocks)
        marker_probit_tree_state = initialize_marker_probit_tree_state(annotation_design, context.startPi)
    end
    block_data = nothing
    GC.gc()

    effect_lengths = effect_lengths_from_blocks(blocks, nCategory)

    if !is_continue
        betaArray = initialize_compact_effect_arrays(nTraits, nCategory, effect_lengths)
        alphaArray = initialize_compact_effect_arrays(nTraits, nCategory, effect_lengths)
        deltaArray = initialize_compact_effect_arrays(nTraits, nCategory, effect_lengths)
    else
        effect_starting_path = starting_value_dir
        delta_starting_path = starting_value_dir * "last_sample_delta/"
        betaArray, alphaArray, deltaArray = initialize_effect_state!(
            effect_starting_path,
            delta_starting_path,
            my_rank,
            my_nsnp,
            nCategory,
            nTraits,
            blocks,
        )
    end

    #posterior mean 
    meanAlpha = initialize_compact_effect_arrays(nTraits, nCategory, effect_lengths)
    
    # mcmc samples for delta -> used to compute PP of SNPs 
    mcmc_Delta = nothing
    if output_mcmc_delta
        nOutput = Int(floor((nIter - burnin) / thin))
        mcmc_Delta = [zeros(my_nsnp * nCategory, nOutput) for _ in 1:nTraits]
    end

    save_category_correlation_outputs = annotation_prior_model == :group_dirichlet

    rank0_mcmc_state = nothing
    rank0_mcmc_state = initialize_rank0_mcmc_state(
        Pi,
        nIter,
        burnin,
        thin,
        nTraits,
        nCategory;
        estimate_pi=estimate_pi,
        estimate_vara=estimate_vara,
        estimate_vare=estimate_vare,
        report_pleiotropic_qtl_effect_matrix=report_pleiotropic_qtl_effect_matrix,
        save_category_correlation_outputs=save_category_correlation_outputs,
    )
    (; mean_pi, mean_pi2, meanB2, meanA2, meanBcor2, meanAcor2, meanA, meanAcor, meanB, meanBcor, meanG, meanG2, meanGcor, meanGcor2, meanGcor_count, meanSSE, meanGtotal, meanGtotal2, meanGcor_total, meanGcor_total2, meanGcor_total_count, mcmcAtruecor_c, mcmcBcor_c, meanR, meanR2) = rank0_mcmc_state

    file_names = nothing
    file_names = prepare_mcmc_output_files(
        analysis_path;
        report_pleiotropic_qtl_effect_matrix=report_pleiotropic_qtl_effect_matrix,
    )
    
    println("---------------- Summary Start --------------")
    println("nIter=$nIter, thin=$thin, seed=$(config.seed), burnin = $burnin")
    println("Julia threads=$(nthreads())")
    println("startPi is: $Pi")
    println("estimate_vare=$estimate_vare,estimate_vara=$estimate_vara")
    println("estimate_pi=$estimate_pi")
    println("analysis_path=$analysis_path")
    println("data_path=$(config.data_path)")
    println("annotation_prior_model=$annotation_prior_model")
    println("report_pleiotropic_qtl_effect_matrix=$report_pleiotropic_qtl_effect_matrix")
    println("output_mcmc_delta=$output_mcmc_delta")
    if estimate_vara
        println("estimate_Gscale=$estimate_Gscale")
        println("starting value of scale_G_vec is: ", scale_G_vec)
    end
    if annotation_prior_model == :group_dirichlet
        println("nCat = $(annotation_metadata.nCat), nCon = $nCon")
    else
        println("effect categories = $nCategory, annotation design features = $(size(marker_probit_tree_state.design_matrix, 2))")
    end
    time_start = now()
    println("Start time: ", time_start)
    println("---------------- Summary End ----------------")

    println("There are $my_nblk LD blocks, and $my_nsnp SNPs in total.")

    nlabel = 4 # number of labels for Pi: (1.0, 1.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0)
    iout = 1 
    iter_after_burnin_thin_index = 1
    last_saved_iter = nIter - (nIter - burnin) % thin
    use_threaded_blocks = nthreads() > 1 && my_nblk > 1
    block_ranges = build_block_ranges(blocks, min(nthreads(), my_nblk))
    block_workspaces = [
        make_nonmpi_block_workspace(blocks, block_range, nCategory, nlabel, nTraits)
        for block_range in block_ranges
    ]

    R_blk = initialize_r_blk_state(
        my_nblk,
        my_rank,
        Rprior;
        is_continue=is_continue,
        starting_value_dir=starting_value_dir,
        fixed_r_dir=estimate_vare ? nothing : config.starting_value_dir,
    )

    estGscale_iter = min(estGscale_iter, nIter)

    # varg Computed as dot(X\alpha, X\alpha)  
    varg_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # varg for different category & trait (saved for hsq computation)
    varg_cov_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # genetic covariance for different category
    # compute total varg for each trait without split into different categories
    totalvarg_blk = [zeros(nTraits, nTraits) for _ in 1:my_nblk]

    # ssq_blk_cat Computed as dot(\alpha, \alpha)
    ssq_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # sum of square for different category & trait (saved for enrichment)

    nLoci_array_vec = [fill(0, nlabel) for c in 1:nCategory] # number of loci (pi) in each category for SNPs in this rank
    # SSE_vec Computed as dot(\beta, \beta) for whole rank
    SSE_vec = [zeros(nTraits, nTraits) for c in 1:nCategory]
    # G_vec Computed as dot(X\alpha, X\alpha) for whole rank
    G_vec = [zeros(nTraits, nTraits) for c in 1:nCategory] # genetic variance for different category for SNPs in this rank

    R_blkmean = zeros(nTraits, nTraits)
    varg_cat_sum = zeros(nTraits, nCategory)
    varg_cov_cat_sum = zeros(nCategory)
    totalvarg = zeros(nTraits, nTraits)
    Atrue_vec = [zeros(nTraits, nTraits) for _ in 1:nCategory]
    nQTL = zeros(Int, nCategory)

    @showprogress "running MCMC ..." for iter = 1:nIter
        
        within_estGscale = estimate_Gscale && iter <= estGscale_iter
        do_thin = iter > burnin && (iter - burnin) % thin == 0

        if do_thin
            iIter = 1.0 / iter_after_burnin_thin_index
        end

        if within_estGscale
            iIter_scaleG = 1.0 / iter
        end

        need_block_totalvarg = estimate_vare || do_thin

        # Reset per-iteration counters that are updated incrementally.
        for c in 1:nCategory
            fill!(nLoci_array_vec[c], 0.0)
        end

        if use_threaded_blocks
            tasks = map(eachindex(block_ranges)) do range_index
                block_range = block_ranges[range_index]
                Threads.@spawn run_nonmpi_block_range!(
                    blocks,
                    block_range,
                    betaArray,
                    alphaArray,
                    deltaArray,
                    meanAlpha,
                    Ainv_vec,
                    Pi,
                    marker_probit_tree_state,
                    R_blk,
                    varg_blk_cat,
                    varg_cov_blk_cat,
                    totalvarg_blk,
                    ssq_blk_cat,
                    block_workspaces[range_index];
                    annotation_prior_model=annotation_prior_model,
                    effect_nCon=effect_nCon,
                    my_nsnp=my_nsnp,
                    nCategory=nCategory,
                    nTraits=nTraits,
                    nlabel=nlabel,
                    do_thin=do_thin,
                    iIter=do_thin ? iIter : 0.0,
                    need_block_totalvarg=need_block_totalvarg,
                    estimate_vare=estimate_vare,
                    df_R=df_R,
                    scale_R=scale_R,
                )
            end
            for task in tasks
                merge_nloci_counts!(nLoci_array_vec, fetch(task))
            end
        else
            merge_nloci_counts!(
                nLoci_array_vec,
                run_nonmpi_block_range!(
                    blocks,
                    1:my_nblk,
                    betaArray,
                    alphaArray,
                    deltaArray,
                    meanAlpha,
                    Ainv_vec,
                    Pi,
                    marker_probit_tree_state,
                    R_blk,
                    varg_blk_cat,
                    varg_cov_blk_cat,
                    totalvarg_blk,
                    ssq_blk_cat,
                    block_workspaces[1];
                    annotation_prior_model=annotation_prior_model,
                    effect_nCon=effect_nCon,
                    my_nsnp=my_nsnp,
                    nCategory=nCategory,
                    nTraits=nTraits,
                    nlabel=nlabel,
                    do_thin=do_thin,
                    iIter=do_thin ? iIter : 0.0,
                    need_block_totalvarg=need_block_totalvarg,
                    estimate_vare=estimate_vare,
                    df_R=df_R,
                    scale_R=scale_R,
                ),
            )
        end

        # get true A matrix by alphaArray
        # use only pleiotropic markers to compute QTL effect variance matrix
        if do_thin && report_pleiotropic_qtl_effect_matrix
            for cat = 1:nCategory
                Atrue_cat = Atrue_vec[cat]
                fill!(Atrue_cat, 0.0)
                nQTL[cat] = 0

                alpha1 = alphaArray[1][cat]
                alpha2 = alphaArray[2][cat]
                for marker in eachindex(alpha1, alpha2)
                    a1 = alpha1[marker]
                    a2 = alpha2[marker]
                    if a1 != 0 && a2 != 0
                        nQTL[cat] += 1
                        Atrue_cat[1, 1] += a1 * a1
                        Atrue_cat[1, 2] += a1 * a2
                        Atrue_cat[2, 2] += a2 * a2
                    end
                end

                Atrue_cat[2, 1] = Atrue_cat[1, 2]
                if nQTL[cat] > 0
                    Atrue_cat ./= nQTL[cat]
                end
                update_matrix_mean!(meanA[cat], Atrue_cat, iIter)
                if estimate_vara 
                    update_matrix_second_moment!(meanA2[cat], Atrue_cat, iIter)
                end
                if save_category_correlation_outputs
                    mcmcAtruecor_c[iter_after_burnin_thin_index, cat] = compute_correlation(Atrue_cat)
                end
                # save the Atrue_cat to file
                open(file_names["marker_effects_variance"], "a") do io
                    writedlm(io, Atrue_cat, ',')
                end
            end
        end

        ########################
        ### Step1. sample Pi ###
        ########################        
        if estimate_pi
            if annotation_prior_model == :group_dirichlet
                tempPi_vec = [zeros(nlabel) for cat = 1:nCategory]
                for cat = 1:nCategory
                    tempPi_vec[cat] = rand(Dirichlet(nLoci_array_vec[cat] .+ 1))
                    if do_thin
                        tempPi2 = tempPi_vec[cat] .^ 2
                        for (iCategori, key) in enumerate(pi_key_order())
                            mean_pi[cat][key] += (tempPi_vec[cat][iCategori] - mean_pi[cat][key]) * iIter
                            mean_pi2[cat][key] += (tempPi2[iCategori] - mean_pi2[cat][key]) * iIter
                        end
                    end
                end

                for cat = 1:nCategory
                    for (iCategori, key) in enumerate(pi_key_order())
                        Pi[cat][key] = tempPi_vec[cat][iCategori] # annotation specific pi
                    end
                end
            else # :marker_probit_tree
                dense_deltaArray = expand_effect_arrays(deltaArray, blocks, my_nsnp, nCategory, nTraits)
                updated_pi = update_marker_probit_tree_priors!(marker_probit_tree_state, dense_deltaArray)
                for (key, value) in updated_pi
                    Pi[1][key] = value
                end
                if do_thin
                    record_marker_probit_tree_coefficient_moments!(marker_probit_tree_state, iIter)
                    for key in pi_key_order()
                        mean_pi[1][key] += (Pi[1][key] - mean_pi[1][key]) * iIter
                        mean_pi2[1][key] += (Pi[1][key]^2 - mean_pi2[1][key]) * iIter
                    end
                end
            end
        end

        ########################################################################
        ### Step2. sample beta effect covariance matrix and Gscale #############
        ########################################################################   
        
        # get SSE_vec (beta'beta) to sample A and to estimate Gscale
        for cat in 1:nCategory
            for traiti in 1:nTraits
                beta_i = betaArray[traiti][cat]
                for traitj in traiti:nTraits
                    beta_j = betaArray[traitj][cat]
                    SSE_vec[cat][traiti, traitj] = dot(beta_i, beta_j)
                    SSE_vec[cat][traitj, traiti] = SSE_vec[cat][traiti, traitj]
                end
            end
        end
        
        # get the running average of SSE to sample Gscale
        if within_estGscale
            for cat = 1:nCategory
                update_matrix_mean!(meanSSE[cat], SSE_vec[cat], iIter_scaleG)
                copyto!(scale_G_vec[cat], meanSSE[cat])
                scale_G_vec[cat] .*= (df_G - nTraits - 1) / nLoci_annot[cat]
            end
            if iter == estGscale_iter
                # save the scale_G_vec
                for cat = 1:nCategory
                    writedlm(analysis_path * "scale_G" * string(cat) * ".txt", scale_G_vec[cat])
                end
            end
        end

        ### sample beta covariance matrix ###
        for cat = 1:nCategory
            if estimate_vara 
                A_vec_sampler_cat = convert(Array, Symmetric(scale_G_vec[cat] + SSE_vec[cat]))
                A_vec[cat] = rand(InverseWishart(df_G + nLoci_annot[cat], A_vec_sampler_cat))
            end
            if do_thin
                # save mean for beta effect variance 
                update_matrix_mean!(meanB[cat], A_vec[cat], iIter)
                if save_category_correlation_outputs
                    mcmcBcor_c[iter_after_burnin_thin_index, cat] = compute_correlation(A_vec[cat])
                end
                if estimate_vara 
                    update_matrix_second_moment!(meanB2[cat], A_vec[cat], iIter)
                end
            end
        end
        
        # update Ainv_vec if needed for sampling beta in next iteration
        if estimate_vara 
            for cat in 1:nCategory
                copyto!(Ainv_vec[cat], inv(A_vec[cat]))
            end
        end

        ########################################################################
        ### Step3. get average residual variance across blocks ##################
        ########################################################################  
        # summing residual variance
        if estimate_vare
            if do_thin
                accumulate_matrices!(R_blkmean, R_blk)
                R_blkmean ./= my_nblk
                update_matrix_mean!(meanR, R_blkmean, iIter)
                update_matrix_second_moment!(meanR2, R_blkmean, iIter)
            end
        end

        ########################################################################################
        ### Step4. get annotation-specific & total genetic variance across blocks ##############
        ########################################################################################
        if do_thin
            accumulate_matrices!(varg_cat_sum, varg_blk_cat) # sum varg_blk_cat to get varg for different category across blocks
            accumulate_vectors!(varg_cov_cat_sum, varg_cov_blk_cat) # sum varg_cov_blk_cat to get genetic covariance for different category across blocks
            accumulate_matrices!(totalvarg, totalvarg_blk) # sum totalvarg_blk to get total genetic variance across all blocks

            for cat = 1:nCategory
                G_vec[cat][1, 1] = varg_cat_sum[1, cat]
                G_vec[cat][2, 2] = varg_cat_sum[2, cat]
                G_vec[cat][1, 2] = G_vec[cat][2, 1] = varg_cov_cat_sum[cat]
            end

            total_gcor = compute_correlation(totalvarg)
            if isfinite(total_gcor)
                meanGcor_total += total_gcor
                meanGcor_total2 += total_gcor^2
                meanGcor_total_count += 1
            end

            for cat = 1:nCategory
                update_matrix_mean!(meanG[cat], G_vec[cat], iIter)
                update_matrix_second_moment!(meanG2[cat], G_vec[cat], iIter)
                if save_category_correlation_outputs
                    gcor = compute_correlation(G_vec[cat])
                    if isfinite(gcor)
                        meanGcor[cat] += gcor
                        meanGcor2[cat] += gcor^2
                        meanGcor_count[cat] += 1
                    end
                end
            end
            update_matrix_mean!(meanGtotal, totalvarg, iIter)
            update_matrix_second_moment!(meanGtotal2, totalvarg, iIter)
        end

        # save MCMC samples & last samples
        if iter > burnin 
            if iter > burnin && (iter - burnin) % thin == 0
                println("iter $iter")
                if output_mcmc_delta
                    for trait = 1:nTraits
                        mcmc_Delta[trait][:, iout] = expand_effect_trait(deltaArray[trait], blocks, my_nsnp, nCategory)
                    end
                end

                for cat = 1:nCategory
                    open(file_names["pi"], "a") do io
                        write_pi_dict(io, Pi[cat])
                    end
                    open(file_names["beta_effects_variance"], "a") do io
                        writedlm(io, A_vec[cat], ',')
                    end
                end
                
                if output_mcmc_delta
                    iout += 1
                end
            end

            if do_thin
                open(file_names["total_genetic_effects_variance"], "a") do io
                    writedlm(io, totalvarg, ',')
                end
                for cat in 1:nCategory
                    open(file_names["genetic_effects_variance"], "a") do io
                        writedlm(io, G_vec[cat], ',')
                    end
                end

                iter_after_burnin_thin_index += 1
            end

            if iter == last_saved_iter
                write_nonmpi_restart_state(
                    analysis_path,
                    my_rank,
                    my_nblk,
                    nCategory,
                    nTraits,
                    expand_effect_arrays(betaArray, blocks, my_nsnp, nCategory, nTraits),
                    R_blk,
                    A_vec,
                    Pi,
                    expand_effect_arrays(deltaArray, blocks, my_nsnp, nCategory, nTraits),
                )
            end
        end

        # GC 
        if iter % 1000 == 0
            GC.gc()
        end
    end # end MCMC iteration loop

    posterior_mean_state = (
        mean_pi=mean_pi,
        mean_pi2=mean_pi2,
        meanB2=meanB2,
        meanA2=meanA2,
        meanBcor2=meanBcor2,
        meanAcor2=meanAcor2,
        meanA=meanA,
        meanAcor=meanAcor,
        meanB=meanB,
        meanBcor=meanBcor,
        meanG=meanG,
        meanG2=meanG2,
        meanGcor=meanGcor,
        meanGcor2=meanGcor2,
        meanGcor_count=meanGcor_count,
        meanGtotal=meanGtotal,
        meanGtotal2=meanGtotal2,
        meanGcor_total=meanGcor_total,
        meanGcor_total2=meanGcor_total2,
        meanGcor_total_count=meanGcor_total_count,
        mcmcAtruecor_c=mcmcAtruecor_c,
        mcmcBcor_c=mcmcBcor_c,
        meanR=meanR,
        meanR2=meanR2,
    )

    save_nonmpi_posterior_mean!(
        analysis_path,
        my_rank,
        nCategory,
        expand_effect_arrays(meanAlpha, blocks, my_nsnp, nCategory, nTraits),
        posterior_mean_state;
        estimate_vara=estimate_vara,
        estimate_vare=estimate_vare,
        estimate_pi=estimate_pi,
        report_pleiotropic_qtl_effect_matrix=report_pleiotropic_qtl_effect_matrix,
        annotation_prior_model=annotation_prior_model,
        marker_probit_tree_state=marker_probit_tree_state,
        annotation_names=annotation_metadata.annotationName,
    )

    if output_mcmc_delta
        writedlm(analysis_path * "mcmc_Delta1.txt", mcmc_Delta[1])
        writedlm(analysis_path * "mcmc_Delta2.txt", mcmc_Delta[2])
    end

    time_end = now()
    time_diff = (time_end - time_start).value / 60000 #milliseconds to min
    println("End time: ", time_end)
    println("Running Time (min): ", time_diff)
    return nothing
end

function run_nonmpi_workflow(config::ConfigTypes.NonMPIConfig)
    Random.seed!(config.seed)
    if nthreads() > 1
        BLAS.set_num_threads(1)
    end
    context = build_nonmpi_run_context(config)
    return @time run_nonmpi_sampler!(context)
end
