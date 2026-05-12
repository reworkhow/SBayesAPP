using DelimitedFiles
using Distributions
using LinearAlgebra
using ProgressMeter
using Random
using Statistics
using Dates
using JLD2


function sample_variance_sumstats(ycorr_array, nobs, df, scale)
    ntraits = length(ycorr_array)
    SSE = zeros(ntraits, ntraits)
    for traiti = 1:ntraits
        ycorri = ycorr_array[traiti]
        for traitj = traiti:ntraits
            ycorrj = ycorr_array[traitj]
            SSE[traiti, traitj] = dot(ycorri, ycorrj)
            SSE[traitj, traiti] = SSE[traiti, traitj]
        end
    end
    R = rand(InverseWishart(df + nobs, convert(Array, Symmetric(scale + SSE))))
    return R
end

function build_nonmpi_sampler_settings(config::ConfigTypes.NonMPIConfig)
    estimate_vara = config.estimate_vara
    estimate_pi = config.estimate_pi
    if config.annotation_prior_model == :marker_probit_tree && !estimate_pi
        @warn "estimate_pi=false is ignored when annotation_prior_model=:marker_probit_tree."
        estimate_pi = true
    end
    is_continue = config.is_continue
    if config.annotation_prior_model == :marker_probit_tree && is_continue
        @warn "is_continue=true is ignored when annotation_prior_model=:marker_probit_tree. A fresh run will be started instead."
        is_continue = false
    end
    return (
        estimate_vare=config.estimate_vare,
        estimate_vara=estimate_vara,
        estimate_pi=estimate_pi,
        estimate_Gscale=config.estimate_Gscale && estimate_vara,
        estGscale_iter=config.estGscale_iter,
        is_continue=is_continue,
    )
end

function build_nonmpi_run_context(config::ConfigTypes.NonMPIConfig)
    annotation_metadata = load_annotation_metadata(config.data_path, config.annot_file; nCon=config.n_con)
    if config.annotation_prior_model == :marker_probit_tree && config.n_con != 0
        @warn "n_con is ignored when annotation_prior_model=:marker_probit_tree; all annotation columns are treated as prior features."
    end
    settings = build_nonmpi_sampler_settings(config)
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

function run_nonmpi_sampler!(context)
    config = context.config
    annotation_metadata = context.annotation_metadata
    settings = context.settings
    Gprior_vec = context.Gprior_vec

    analysis_path = config.analysis_path
    nIter = config.nIter
    outFreq = config.out_freq
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
    if config.nrank == 1
        my_rank = 0
    end

    block_data = load_nonmpi_block_data(config.data_path, config.annot_dict)
    my_TransformedY_dict = block_data.transformed_y_dict
    my_blkSNPsIndex_dict = block_data.blkSNPsIndex_dict
    my_blkID = block_data.blkID
    my_nGWAS_dict = block_data.nGWAS_dict
    my_nblk = block_data.nblk
    my_nsnp = block_data.nsnp
    raw_anno_matrix_dict = block_data.anno_matrix_dict

    if is_continue
        burnin = 0
    else 
        burnin = floor(Int, nIter * 0.4)
    end

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

    if my_rank == 0
        writedlm(analysis_path * "annotationName.txt", annotation_metadata.annotationName)
    end
    marker_probit_tree_state = nothing
    if annotation_prior_model == :group_dirichlet
        xpx_dict, xArray_dict, my_anno_matrix_dict = prepare_block_state!(
            my_blkID,
            my_blkSNPsIndex_dict,
            block_data.transformed_x_dict,
            raw_anno_matrix_dict,
            effect_nCon,
            annotation_metadata.annotationType,
        )
    else
        xpx_dict, xArray_dict, my_anno_matrix_dict, annotation_design = prepare_marker_probit_tree_block_state!(
            my_blkID,
            my_blkSNPsIndex_dict,
            block_data.transformed_x_dict,
            raw_anno_matrix_dict,
        )
        marker_probit_tree_state = initialize_marker_probit_tree_state(annotation_design, context.startPi)
    end
    block_data = nothing
    GC.gc()

    #output
    β = zeros(nTraits)
    newα = zeros(nTraits)  #α=Dβ
    oldα = zeros(nTraits)
    w = zeros(nTraits)
    δ = zeros(nTraits)

    if !is_continue
        betaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits] #-> ordered by annotation groups 
        alphaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
        deltaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
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
            my_blkID,
            my_blkSNPsIndex_dict,
            xArray_dict,
            my_TransformedY_dict,
        )
    end

    #posterior mean 
    meanAlpha = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
    
    # mcmc samples for delta -> used to compute PP of SNPs 
    mcmc_Delta = nothing
    if output_mcmc_delta
        nOutput = Int(floor((nIter - burnin) / outFreq))
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
    (; mean_pi, mean_pi2, meanB2, meanA2, meanBcor2, meanAcor2, meanA, meanAcor, meanB, meanBcor, meanG, meanG2, meanGcor, meanGcor2, meanSSE, meanGtotal, meanGtotal2, mcmcAtruecor_c, mcmcBcor_c, mcmcGcov_c, mcmcGcor_c, mcmcGcov_total, mcmcGcor_total, meanR, meanR2) = rank0_mcmc_state

    file_names = nothing
    if my_rank == 0
        file_names = prepare_mcmc_output_files(
            analysis_path;
            report_pleiotropic_qtl_effect_matrix=report_pleiotropic_qtl_effect_matrix,
        )
    end

    if my_rank == 0
        println("---------------- Summary Start --------------")
        println("nIter=$nIter, outFreq=$outFreq, seed=$(config.seed), burnin = $burnin")
        println("startPi is: $Pi")
        println("Number of ranks: ", config.nrank)
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
        println("thin = $thin")
        time_start = now()
        println("Start time: ", time_start)
        println("---------------- Summary End ----------------")
    end
    println("In rank$my_rank, there are $my_nblk LD blocks, and $my_nsnp SNPs in total.")

    nlabel = 4 # number of labels for Pi: (1.0, 1.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0)
    iout = 1 
    iter_after_burnin_thin_index = 1
    last_saved_iter = nIter - (nIter - burnin) % outFreq

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

        #for loop for LD blocks
        for b in 1:my_nblk
            blk = my_blkID[b] #block ID
            xpx_vec = xpx_dict[blk]
            xArray_vec = xArray_dict[blk]
            wArray = my_TransformedY_dict[blk]
            annotationMatb = my_anno_matrix_dict[blk]  #annotation boolean data for current blk, nMarkerb-by-C
            SNPIndexb = my_blkSNPsIndex_dict[blk]
            nEigenb, nMarkerb = size(xArray_vec[end]) # q & nsnpb
            nInd = my_nGWAS_dict[blk] # n

            # # initialize xArrayc/xpxc 
            xArrayc = xArray_vec[1] # nEigenb x nMarkerb matrix
            xpxc = xpx_vec[1]

            Rinv = inv(R_blk[b])

            MarkerOrder = shuffle(1:nMarkerb)
            for marker = MarkerOrder
                true_marker_num = SNPIndexb[marker] #marker position across my_nsnp
                annoindexm = annotationMatb[marker, :] # categories for current marker

                for cat = 1:nCategory
                    if effect_nCon != 0
                        if cat <= effect_nCon + 1 # continuous group + categorical group (after effect_nCon+1, xArrayc/xpxc will not change)
                            xArrayc = xArray_vec[cat]
                            xpxc = xpx_vec[cat]
                        end
                    end
                    Ginv = Ainv_vec[cat]

                    if annoindexm[cat] # skip sampling if the marker is not in the category
                        markerIndex = (cat - 1) * my_nsnp + true_marker_num # exact position in betaArray; betaArray[trait1]: nMarker*nCategory-by-1
                        x = xArrayc[:, marker]
                        for trait = 1:nTraits
                            β[trait] = betaArray[trait][markerIndex]
                            oldα[trait] = newα[trait] = alphaArray[trait][markerIndex]
                            δ[trait] = deltaArray[trait][markerIndex]
                            w[trait] = dot(x, wArray[trait]) + xpxc[marker] * oldα[trait] #w=xj'(ycorr+xj*αj) scaler
                        end

                        for k = 1:nTraits
                            Ginv11 = Ginv[k, k]
                            nok = deleteat!(collect(1:nTraits), k)
                            Ginv12 = Ginv[k, nok] #this is not row vector!!, so this is Ginv21
                            C11 = Ginv11 + Rinv[k, k] * xpxc[marker]
                            C12 = Ginv12 + xpxc[marker] * Diagonal(δ[nok]) * Rinv[k, nok] #C21
                            #when δj=0
                            invLhs0 = 1 / Ginv11
                            rhs0 = -dot(Ginv12, β[nok])
                            gHat0 = rhs0 * invLhs0
                            #when δj=1
                            invLhs1 = 1 / C11
                            rhs1 = w' * Rinv[:, k] - C12'β[nok]  #here the w' is in paper: w'm=xj'(ycorr+xj*βj)
                            gHat1 = rhs1 * invLhs1

                            d0 = deepcopy(δ)
                            d1 = deepcopy(δ)
                            d0[k] = 0.0
                            d1[k] = 1.0

                            #sample δj
                            logPrior0 = log_marker_state_prior(annotation_prior_model, Pi, marker_probit_tree_state, cat, true_marker_num, d0)
                            logPrior1 = log_marker_state_prior(annotation_prior_model, Pi, marker_probit_tree_state, cat, true_marker_num, d1)
                            logDelta0 = -0.5 * (log(Ginv11) - gHat0^2 * Ginv11) + logPrior0 #logPi
                            logDelta1 = -0.5 * (log(C11) - gHat1^2 * C11) + logPrior1 #logPiComp
                            probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1))

                            #sample marker effects
                            if (rand() < probDelta1) #δj=1
                                δ[k] = 1
                                β[k] = newα[k] = gHat1 + randn() * sqrt(invLhs1)
                                wArray[k] = wArray[k] + x * (oldα[k] - newα[k])
                            else
                                β[k] = gHat0 + randn() * sqrt(invLhs0)
                                δ[k] = 0
                                newα[k] = 0
                                if oldα[k] != 0
                                    wArray[k] = wArray[k] + x * oldα[k] #newα[k]=0
                                end
                            end
                        end

                        # add to nLoci_array_vec based on δ
                        state_key = pi_key(δ)
                        for (pi_index, key) in enumerate(pi_key_order())
                            if state_key == key
                                nLoci_array_vec[cat][pi_index] += 1
                                break
                            end
                        end

                        for trait = 1:nTraits
                            betaArray[trait][markerIndex] = β[trait]
                            deltaArray[trait][markerIndex] = δ[trait]
                            alphaArray[trait][markerIndex] = newα[trait]
                            if do_thin
                                meanAlpha[trait][markerIndex] += (newα[trait] - meanAlpha[trait][markerIndex]) * iIter
                            end
                        end
                    end # end if loop  
                end # end annotation loop
            end  # end marker loop

            my_TransformedY_dict[blk] = wArray
            
            if need_block_totalvarg
                # compute genetic variance for each category and total genetic variance for current block based on alphaArray
                XAb = xArray_vec[1]
                what_array_total = [zeros(size(XAb, 1)) for _ in 1:nTraits] # used to compute total genetic variance

                for cat in 1:nCategory
                    if effect_nCon != 0 && cat <= effect_nCon + 1 # continuous group + categorical group (after effect_nCon+1, xArrayc/xpxc will not change)
                        XAb = xArray_vec[cat]
                    end
                    alphaArray_c = [alphaArray[i][(cat-1)*my_nsnp.+SNPIndexb] for i in 1:nTraits]
                    what_array_c = [XAb * alphaArray_c[traiti] for traiti in 1:nTraits]
                    for traiti = 1:nTraits
                        if do_thin
                            varg_blk_cat[b][traiti, cat] = dot(what_array_c[traiti], what_array_c[traiti]) # genetic variance for different category & trait in bth block
                        end
                        if estimate_vare
                            ssq_blk_cat[b][traiti, cat] = dot(alphaArray_c[traiti], alphaArray_c[traiti])
                        end
                        what_array_total[traiti] += what_array_c[traiti]  # add up whati across category to compute total genetic variance
                    end
                    if do_thin
                        varg_cov_blk_cat[b][cat] = dot(what_array_c[1], what_array_c[2]) # genetic covariance for different category in bth block
                    end
                end

                # total genetic variance
                for traiti in 1:nTraits
                    for traitj in traiti:nTraits
                        totalvarg_blk[b][traiti, traitj] = dot(what_array_total[traiti], what_array_total[traitj])
                        totalvarg_blk[b][traitj, traiti] = totalvarg_blk[b][traiti, traitj]
                    end
                end
            end
           
            if estimate_vare
                # sample bivariate residual variance R 
                sampled_R = sample_variance_sumstats(wArray, nEigenb, df_R, scale_R)
                R_blk[b] = sampled_R
                Rcor = compute_correlation(sampled_R)
                for traiti = 1:nTraits
                    thres = sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]
                    if (thres > 1.1)
                        R_blk[b][traiti, traiti] = sampled_R[traiti, traiti]
                    else
                        R_blk[b][traiti, traiti] = 1.0 / nInd[traiti]
                    end
                end
                # tune covariance in R_blk
                Rcov = Rcor * sqrt(R_blk[b][1, 1] * R_blk[b][2, 2])
                R_blk[b][1, 2] = R_blk[b][2,1] = Rcov     
            end 
        end # end block loop 

        # get true A matrix by alphaArray
        # use only pleiotropic markers to compute QTL effect variance matrix
        if do_thin && report_pleiotropic_qtl_effect_matrix
            Atrue_vec = [zeros(nTraits, nTraits) for c in 1:nCategory]
            nQTL = zeros(nCategory)
            for cat = 1:nCategory
                alpha_array = [alphaArray[i][((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)] for i in 1:nTraits]
                pleio_marker = (alpha_array[1] .!= 0) .& (alpha_array[2] .!= 0)
                nQTL[cat] = sum(pleio_marker)
                for traiti in 1:nTraits
                    for traitj in traiti:nTraits
                        Atrue_vec[cat][traiti, traitj] = dot(alpha_array[traiti][pleio_marker], alpha_array[traitj][pleio_marker])
                        Atrue_vec[cat][traitj, traiti] = Atrue_vec[cat][traiti, traitj]
                    end
                end
            end
            if my_rank == 0
                for cat = 1:nCategory
                    #annotation specific A
                    if nQTL[cat] == 0
                        Atrue_cat = zeros(nTraits, nTraits)
                    else
                        Atrue_cat = Atrue_vec[cat] / nQTL[cat]
                    end
                    meanA[cat] += (Atrue_cat - meanA[cat]) * iIter
                    if estimate_vara 
                        meanA2[cat] += (Atrue_cat .^ 2 - meanA2[cat]) * iIter
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
                        Pi[cat][key] = tempPi_vec[cat][iCategori] #annotation specific pi
                    end
                end
            else # :marker_probit_tree
                updated_pi = update_marker_probit_tree_priors!(marker_probit_tree_state, deltaArray)
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
                beta_i = betaArray[traiti][((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)]
                for traitj in traiti:nTraits
                    beta_j = betaArray[traitj][((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)]
                    SSE_vec[cat][traiti, traitj] = dot(beta_i, beta_j)
                    SSE_vec[cat][traitj, traiti] = SSE_vec[cat][traiti, traitj]
                end
            end
        end
        
        # get the running average of SSE to sample Gscale
        if within_estGscale
            if my_rank == 0
                for cat = 1:nCategory
                    meanSSE[cat] += (SSE_vec[cat] - meanSSE[cat]) * iIter_scaleG
                end
            end

            Gprior_vec = deepcopy(meanSSE)
            for cat = 1:nCategory
                Gprior_vec[cat] = meanSSE[cat] / nLoci_annot[cat]
                scale_G_vec[cat] = Gprior_vec[cat] * (df_G - nTraits - 1)
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
                A_vec_sampler = convert(Array, Symmetric(scale_G_vec[cat] + SSE_vec[cat]))
                A_vec[cat] = rand(InverseWishart(df_G + nLoci_annot[cat], A_vec_sampler))  
            end
            if do_thin
                # save mean for beta effect variance 
                meanB[cat] += (A_vec[cat] - meanB[cat]) * iIter
                if save_category_correlation_outputs
                    mcmcBcor_c[iter_after_burnin_thin_index, cat] = compute_correlation(A_vec[cat])
                end
                if estimate_vara 
                    meanB2[cat] += (A_vec[cat] .^ 2 - meanB2[cat]) * iIter
                end
            end
        end
        
        # update Ainv_vec if needed for sampling beta in next iteration
        if estimate_vara 
            Ainv_vec[:] = [inv(A_vec[cat]) for cat in 1:nCategory]
        end
        # gc after sampling A
        GC.gc()

        ########################################################################
        ### Step3. get average residual variance across blocks ##################
        ########################################################################  
        # summing residual variance
        if estimate_vare
            if do_thin
                R_blk_sum = sum(R_blk)
                R_blkmean = R_blk_sum / my_nblk
                R2 = (R_blkmean) .^ 2
                meanR += (R_blkmean - meanR) * iIter
                meanR2 += (R2 - meanR2) * iIter
            end
        end

        ########################################################################################
        ### Step4. get annotation-specific & total genetic variance across blocks ##############
        ########################################################################################
        if do_thin
            varg_cat = sum(varg_blk_cat) # sum varg_blk_cat to get varg for different category across blocks
            varg_cov_cat = sum(varg_cov_blk_cat) # sum varg_cov_blk_cat to get genetic covariance for different category across blocks
            totalvarg = sum(totalvarg_blk) # sum totalvarg_blk to get total genetic variance across all blocks

            for cat = 1:nCategory
                G_vec[cat][1, 1] = varg_cat[1, cat]
                G_vec[cat][2, 2] = varg_cat[2, cat]
                G_vec[cat][1, 2] = G_vec[cat][2, 1] = varg_cov_cat[cat]
            end

            mcmcGcor_total[iter_after_burnin_thin_index] = compute_correlation(totalvarg)
            mcmcGcov_total[iter_after_burnin_thin_index] = totalvarg[1, 2]

            for cat = 1:nCategory
                meanG[cat] += (G_vec[cat] - meanG[cat]) * iIter
                meanG2[cat] += (G_vec[cat] .^ 2 - meanG2[cat]) * iIter
                if save_category_correlation_outputs
                    mcmcGcor_c[iter_after_burnin_thin_index, cat] = compute_correlation(G_vec[cat])
                    mcmcGcov_c[iter_after_burnin_thin_index, cat] = G_vec[cat][1, 2]
                end
            end
            meanGtotal += (totalvarg - meanGtotal) * iIter
            meanGtotal2 += (totalvarg .^ 2 - meanGtotal2) * iIter
        end

        # save MCMC samples & last samples
        if iter > burnin 
            if iter > burnin && (iter - burnin) % outFreq == 0
                println("iter $iter")
                if output_mcmc_delta
                    for trait = 1:nTraits
                        mcmc_Delta[trait][:, iout] = deltaArray[trait]
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
                if my_rank == 0
                    open(file_names["total_genetic_effects_variance"], "a") do io
                        writedlm(io, totalvarg, ',')
                    end
                    for cat in 1:nCategory
                        open(file_names["genetic_effects_variance"], "a") do io
                            writedlm(io, G_vec[cat], ',')
                        end
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
                    betaArray,
                    R_blk,
                    A_vec,
                    Pi,
                    deltaArray,
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
        meanGtotal=meanGtotal,
        meanGtotal2=meanGtotal2,
        mcmcAtruecor_c=mcmcAtruecor_c,
        mcmcBcor_c=mcmcBcor_c,
        mcmcGcov_c=mcmcGcov_c,
        mcmcGcor_c=mcmcGcor_c,
        mcmcGcov_total=mcmcGcov_total,
        mcmcGcor_total=mcmcGcor_total,
        meanR=meanR,
        meanR2=meanR2,
    )

    save_nonmpi_posterior_mean!(
        analysis_path,
        my_rank,
        nCategory,
        meanAlpha,
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
        writedlm(analysis_path * "mcmc_Delta1.rank$my_rank.txt", mcmc_Delta[1])
        writedlm(analysis_path * "mcmc_Delta2.rank$my_rank.txt", mcmc_Delta[2])
    end

    time_end = now()
    time_diff = (time_end - time_start).value / 60000 #milliseconds to min
    println("End time: ", time_end)
    println("Running Time (min): ", time_diff)
    return nothing
end

function run_nonmpi_workflow(config::ConfigTypes.NonMPIConfig)
    context = build_nonmpi_run_context(config)
    return @time run_nonmpi_sampler!(context)
end
