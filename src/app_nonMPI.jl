using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP

    ############################################################################
    #read data in current rank
    ############################################################################
    if nrank == 1
        my_rank = 0
    end

    my_TransformedX_dict = JLD2.load(data_path * "TransformedX_dict.jld2")["my_TransformedX_dict"] # dictionary of Q for each block; blk => Q
    my_TransformedY_dict = JLD2.load(data_path * "TransformedY_dict.jld2")["my_TransformedY_dict"] # dictionary of w for each block; blk => [w1,w2]
    my_blkSNPsIndex_dict = JLD2.load(data_path * "blkSNPsIndex_dict.jld2")["my_blkSNPsIndex_dict"] # dictionary of SNPs index for each block; blk => [1,2,...,nsnp_blk]
    my_blkID = Int.(vec(readdlm(data_path * "blkIDs.txt", ',')))  # block IDs to be analyzed 
    my_nGWAS_dict = JLD2.load(data_path * "nGWAS_dict.jld2")["my_nGWAS_dict"] # dictionary of average nInd for each block; blk => [nInd1, nInd2]
    my_nblk = length(my_blkID) # number of blocks analyzed 
    my_nsnp = sum(map(length, values(my_blkSNPsIndex_dict))) # total number of SNPs to be analyzed  
    my_anno_matrix_dict = JLD2.load(data_path * "$annot_dict.jld2")["my_anno_matrix_dict"] # dictionary of annotation matrix for each block; blk => A 

    if is_continue
        burnin = 0
    else 
        burnin = floor(Int, nIter * 0.4)
    end

    # hyper-parameters for A (marker effect variance) and R
    nCategory = nCat + nCon
    # marker effect variance matrix
    A_vec = [zeros(2, 2) for c in 1:nCategory]
    if is_continue
        A_vec_starting_path = starting_value_dir * "beta_effect_var_matrices_last_sample/"
        for c in 1:nCategory
            A_vec[c] = readdlm(A_vec_starting_path * "beta_effect_matrix_$c.txt", ',')
        end
    else # starting values 
        A_vec = Gprior_vec
    end
    Ainv_vec = [inv(A_vec[c]) for c in 1:nCategory]

    if is_continue
        Pi_starting_path = starting_value_dir * "pi_last_sample/"
        Pi = [Dict{NTuple{2,Float64},Float64}() for c in 1:nCategory]
        for c in 1:nCategory
            Pi[c] = read_to_dict(Pi_starting_path * "pi_$c.txt")
        end
    else
        Pi = [deepcopy(startPi) for c in 1:nCategory]
    end

    # starting value of residual variance matrix R
    Rprior = [1. 0. ; 0. 1.]
    nTraits = 2
    if estimate_vare == true
        df_R = 4 + nTraits
        scale_R = Rprior * (df_R - nTraits - 1)
    end

    if my_rank == 0
        if estimate_vara == true
            df_G = 4 + nTraits

            if is_continue && estimate_Gscale
                scale_G_vec = [zeros(2, 2) for c in 1:nCategory]
                for c in 1:nCategory
                    scale_G_vec[c] = readdlm(gscale_value_dir * "scale_G$c.txt")
                end
                estimate_Gscale = false # estimate_Gscale = false for continuous analysis
                println("scale_G_vec is fixed as saved in $gscale_value_dir")
            else
                scale_G_vec = [Gprior_vec[c] * (df_G - nTraits - 1) for c in 1:nCategory]
                println("scale_G_vec is computed by ST h2.")
            end
        end
    end

    if my_rank == 0
        writedlm(analysis_path * "annotationName.txt", annotationName)
    end
    xpx_dict, xArray_dict, my_anno_matrix_dict = prepare_block_state!(
        my_blkID,
        my_blkSNPsIndex_dict,
        my_TransformedX_dict,
        my_anno_matrix_dict,
        nCon,
        annotationType,
    )

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
    nOutput = Int(floor((nIter-burnin) / outFreq))
    mcmc_Delta = [zeros(my_nsnp * nCategory, nOutput) for t in 1:nTraits]

    if my_rank == 0
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
        )
        nsample4mean = rank0_mcmc_state.nsample4mean
        mean_pi = rank0_mcmc_state.mean_pi
        mean_pi2 = rank0_mcmc_state.mean_pi2
        meanB2 = rank0_mcmc_state.meanB2
        meanA2 = rank0_mcmc_state.meanA2
        meanBcor2 = rank0_mcmc_state.meanBcor2
        meanAcor2 = rank0_mcmc_state.meanAcor2
        meanA = rank0_mcmc_state.meanA
        meanAcor = rank0_mcmc_state.meanAcor
        meanB = rank0_mcmc_state.meanB
        meanBcor = rank0_mcmc_state.meanBcor
        meanG = rank0_mcmc_state.meanG
        meanG2 = rank0_mcmc_state.meanG2
        meanGcor = rank0_mcmc_state.meanGcor
        meanGcor2 = rank0_mcmc_state.meanGcor2
        meanSSE = rank0_mcmc_state.meanSSE
        meanGtotal = rank0_mcmc_state.meanGtotal
        meanGtotal2 = rank0_mcmc_state.meanGtotal2
        mcmcAtruecor_c = rank0_mcmc_state.mcmcAtruecor_c
        mcmcBcor_c = rank0_mcmc_state.mcmcBcor_c
        mcmcGcov_c = rank0_mcmc_state.mcmcGcov_c
        mcmcGcor_c = rank0_mcmc_state.mcmcGcor_c
        mcmcGcov_total = rank0_mcmc_state.mcmcGcov_total
        mcmcGcor_total = rank0_mcmc_state.mcmcGcor_total
        meanR = rank0_mcmc_state.meanR
        meanR2 = rank0_mcmc_state.meanR2
    end

    file_names = nothing
    if my_rank == 0
        file_names = prepare_mcmc_output_files(analysis_path)
    end

    if my_rank == 0
        println("---------------- Summary Start --------------")
        println("nIter=$nIter, outFreq=$outFreq, seed=$seed, burnin = $burnin")
        println("startPi is: $Pi")
        println("Number of ranks: ", nrank)
        println("estimate_vare=$estimate_vare,estimate_vara=$estimate_vara")
        println("estimate_pi=$estimate_pi")
        println("estimate_Gscale=$estimate_Gscale")
        println("analysis_path=$analysis_path")
        println("data_path=$data_path")
        if estimate_vara
            println("scale_G_vec is: ", scale_G_vec)
        end
        println("nCat = $nCat, nCon = $nCon")
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
    )

    estGscale_iter = min(estGscale_iter, nIter)

    # varg Computed as dot(X\alpha, X\alpha)  
    varg_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # varg for different category & trait (saved for hsq computation)
    varg_cov_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # genetic covariance for different category
    # compute total varg for each trait without split into different categories
    totalvarg_blk = [zeros(nTraits, nTraits) for _ in 1:my_nblk]

    # ssq_blk_cat Computed as dot(\alpha, \alpha)
    ssq_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # sum of square for different category & trait (saved for enrichment)
    ssq_cov_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # sum of square for different category (saved for enrichment)

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

            R = deepcopy(R_blk[b])
            for traiti = 1:nTraits
                for traitj = traiti:nTraits
                    R[traiti, traitj] = R[traiti, traitj] / sqrt(nInd[traiti] * nInd[traitj])
                    R[traitj, traiti] = R[traiti, traitj]
                end
            end
            Rinv = inv(R)

            MarkerOrder = shuffle(1:nMarkerb)
            for marker = MarkerOrder
                true_marker_num = SNPIndexb[marker] #marker position across my_nsnp
                annoindexm = annotationMatb[marker, :] # categories for current marker

                for cat = 1:nCategory
                    if nCon != 0
                        if cat <= nCon + 1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                            xArrayc = xArray_vec[cat]
                            xpxc = xpx_vec[cat]
                        end
                    end
                    Ginv = Ainv_vec[cat]
                    BigPi = Pi[cat]

                    if annoindexm[cat] # skip sampling if the marker is not in the category
                        markerIndex = (cat - 1) * my_nsnp + true_marker_num # exact position in betaArray; betaArray[trait1]: nMarker*nCategory-by-1
                        x = xArrayc[:, marker]
                        for trait = 1:nTraits
                            β[trait] = betaArray[trait][markerIndex]
                            oldα[trait] = newα[trait] = alphaArray[trait][markerIndex]
                            δ[trait] = deltaArray[trait][markerIndex]
                            w[trait] = dot(x, wArray[trait]) + xpxc[marker] * oldα[trait] #w=xj'(ycorr+xj*αj) scaler
                        end

                        TraitOrder = shuffle(1:nTraits)
                        for k = TraitOrder
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
                            logDelta0 = -0.5 * (log(Ginv11) - gHat0^2 * Ginv11) + log(BigPi[SBayesAPP.pi_key(d0)]) #logPi
                            logDelta1 = -0.5 * (log(C11) - gHat1^2 * C11) + log(BigPi[SBayesAPP.pi_key(d1)]) #logPiComp
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
                        state_key = SBayesAPP.pi_key(δ)
                        for (pi_index, key) in enumerate(SBayesAPP.pi_key_order())
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
                # compute genetic variance and heritability 
                XAb = xArray_vec[1]
                what_array_total = [zeros(size(XAb, 1)) for _ in 1:nTraits] # used to compute total genetic variance

                for cat in 1:nCategory
                    if nCon != 0 && cat <= nCon + 1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                        XAb = xArray_vec[cat]
                    end
                    alphaArray_c = [alphaArray[i][(cat-1)*my_nsnp.+SNPIndexb] for i in 1:nTraits]
                    what_array_c = [XAb * alphaArray_c[traiti] for traiti in 1:nTraits]
                    for traiti = 1:nTraits
                        if do_thin
                            varg_blk_cat[b][traiti, cat] = dot(what_array_c[traiti], what_array_c[traiti]) # genetic variance for different category & trait in bth block
                        end
                        if estimate_vare == true
                            ssq_blk_cat[b][traiti, cat] = dot(alphaArray_c[traiti], alphaArray_c[traiti])
                        end
                        what_array_total[traiti] += what_array_c[traiti]  # add up whati across category to compute total genetic variance
                    end
                    if do_thin
                        varg_cov_blk_cat[b][cat] = dot(what_array_c[1], what_array_c[2]) # genetic covariance for different category in bth block
                        ssq_cov_blk_cat[b][cat] = dot(alphaArray_c[1], alphaArray_c[2])
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
           
            if estimate_vare == true
                # sample bivariate residual variance R 
                sampled_R = sample_variance_sumstats(wArray, nEigenb, df_R, scale_R, nInd)
                R_blk[b] = sampled_R
                Rcor = compute_correlation(sampled_R)
                for traiti = 1:nTraits
                    thres = sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]
                    if (thres > 1.1)
                        R_blk[b][traiti, traiti] = sampled_R[traiti, traiti]
                    else
                        R_blk[b][traiti, traiti] = 1.0
                    end
                end
                # tune covariance in R_blk
                Rcov = Rcor * sqrt(R_blk[b][1, 1] * R_blk[b][2, 2])
                R_blk[b][1, 2] = R_blk[b][2,1] = Rcov     
            end 
        end # end block loop 

        # summing residual variance
        if estimate_vare == true
            R_blk_sum = sum(R_blk)
        end

        # get true A matrix by alphaArray
        # use only pleiotropic markers to compute QTL effect variance matrix
        if do_thin 
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
                    if estimate_vara == true
                        meanA2[cat] += (Atrue_cat .^ 2 - meanA2[cat]) * iIter
                    end
                    # correlation values
                    mcmcAtruecor_c[iter_after_burnin_thin_index, cat] = compute_correlation(Atrue_cat)
                    # save the Atrue_cat to file
                    open(file_names["marker_effects_variance"], "a") do io
                        writedlm(io, Atrue_cat, ',')
                    end
                end
            end
        end

        ########################
        ### Step2. sample Pi ###
        ########################        
        if estimate_pi == true
            tempPi_vec = [zeros(nlabel) for cat = 1:nCategory]
            for cat = 1:nCategory
                tempPi_vec[cat] = rand(Dirichlet(nLoci_array_vec[cat] .+ 1))
                if do_thin
                    tempPi2 = tempPi_vec[cat] .^ 2
                    for (iCategori, key) in enumerate(SBayesAPP.pi_key_order())
                        mean_pi[cat][key] += (tempPi_vec[cat][iCategori] - mean_pi[cat][key]) * iIter
                        mean_pi2[cat][key] += (tempPi2[iCategori] - mean_pi2[cat][key]) * iIter
                    end
                end   
            end
            
            ############################################################
            # reformat the tempPi_vec into the Pi dictionary
            ############################################################
            for cat = 1:nCategory
                for (iCategori, key) in enumerate(SBayesAPP.pi_key_order())
                    Pi[cat][key] = tempPi_vec[cat][iCategori] #annotation specific pi
                end
            end
        end

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
                # correlation values
                mcmcGcor_c[iter_after_burnin_thin_index, cat] = compute_correlation(G_vec[cat])
                mcmcGcov_c[iter_after_burnin_thin_index, cat] = G_vec[cat][1, 2]
                meanG[cat] += (G_vec[cat] - meanG[cat]) * iIter
                meanG2[cat] += (G_vec[cat] .^ 2 - meanG2[cat]) * iIter
            end
            meanGtotal += (totalvarg - meanGtotal) * iIter
            meanGtotal2 += (totalvarg .^ 2 - meanGtotal2) * iIter
        end

        # get SSE_vec to sample A and to estGscale
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
        
        # summing genetic variance
        if within_estGscale
            if my_rank == 0
                for cat = 1:nCategory
                    meanSSE[cat] += (SSE_vec[cat] - meanSSE[cat]) * iIter_scaleG
                end
            end
        end

        ########################
        ### sample scale_G ###
        ########################
        # use empirical SSE to estimate scale_G_vec
        if within_estGscale
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

        ### Step3. sample variance ###
        for cat = 1:nCategory
            if estimate_vara == true
                A_vec_sampler = convert(Array, Symmetric(scale_G_vec[cat] + SSE_vec[cat]))
                A_vec[cat] = rand(InverseWishart(df_G + nLoci_annot[cat], A_vec_sampler))  
            end
            if do_thin
                # save mean for beta effect variance 
                meanB[cat] += (A_vec[cat] - meanB[cat]) * iIter
                # correlation values
                mcmcBcor_c[iter_after_burnin_thin_index, cat] = compute_correlation(A_vec[cat])
                if estimate_vara == true
                    meanB2[cat] += (A_vec[cat] .^ 2 - meanB2[cat]) * iIter
                end
            end
        end
        
        if estimate_vara == true 
            Ainv_vec[:] = [inv(A_vec[cat]) for cat in 1:nCategory]
        end
        
        # sampling residual variance
        if estimate_vare
            R_blkmean = R_blk_sum / my_nblk
            if do_thin
                R2 = (R_blkmean) .^ 2
                meanR += (R_blkmean - meanR) * iIter
                meanR2 += (R2 - meanR2) * iIter
            end
        
            # update block-wize R by R_blkmean for next iteration
            for b in 1:my_nblk
                R_blk[b] = R_blkmean
            end
        end


        # save MCMC samples & last samples
        if iter > burnin 
            if iter - burnin > 0 && (iter - burnin) % outFreq == 0
                println("iter $iter")
                for trait = 1:nTraits
                    mcmc_Delta[trait][:, iout] = deltaArray[trait]
                end

                for cat = 1:nCategory
                    open(file_names["pi"], "a") do io
                        SBayesAPP.write_pi_dict(io, Pi[cat])
                    end
                    open(file_names["beta_effects_variance"], "a") do io
                        writedlm(io, A_vec[cat], ',')
                    end
                end
                
                iout += 1
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
                for t in 1:nTraits
                    writedlm(analysis_path * "last_mcmc_betaArray$t.rank$my_rank.txt", betaArray[t])
                end
                # Save R_blk
                mkpath(analysis_path * "last_sample_R_blk/")  # Create folder if not exist
                for b in 1:my_nblk
                    writedlm(analysis_path * "last_sample_R_blk/R_blk_b$(b)_rank$(my_rank).txt", R_blk[b])
                end
                # Save beta_effect variance matrices (A_vec)
                mkpath(analysis_path * "beta_effect_var_matrices_last_sample/")
                for c in 1:nCategory
                     writedlm(analysis_path * "beta_effect_var_matrices_last_sample/beta_effect_matrix_$(c).txt", A_vec[c], ',')
                end
                # Save Pi from last sample
                mkpath(analysis_path * "pi_last_sample/")
                for c in 1:nCategory
                    SBayesAPP.write_pi_dict(analysis_path * "pi_last_sample/pi_$(c).txt", Pi[c])
                end
                ### Save delta (last column of mcmc_Delta)
                mkpath(analysis_path * "last_sample_delta/")
                for t in 1:nTraits
                    last_delta = mcmc_Delta[t][:, end]  # extract last column
                    open(analysis_path * "last_sample_delta/last_sample_delta$(t)_rank$(my_rank).txt", "w") do io
                        writedlm(io, last_delta)
                    end
                end
            end
        end

        # GC 
        if iter % 1000 == 0
            GC.gc()
        end
    end # end MCMC iteration loop
    
    # save posterior mean
    for cat in 1:nCategory
        # marker effect variance 
        writedlm(analysis_path * "estA" * string(cat) * ".txt", meanA[cat])
        meanAcor[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])])
        # beta effect variance
        writedlm(analysis_path * "estB" * string(cat) * ".txt", meanB[cat])
        meanBcor[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])])
        if estimate_vara == true 
            writedlm(analysis_path * "estB_std" * string(cat) * ".txt", sqrt.((meanB2[cat] .- (meanB[cat] .^ 2))))
            writedlm(analysis_path * "estA_std" * string(cat) * ".txt", sqrt.((meanA2[cat] .- (meanA[cat] .^ 2))))
            meanBcor2[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])] .^ 2)
            meanAcor2[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])].^2)
        end
    end

    # save mcmcGcov_c and mcmcGcov_total
    writedlm(analysis_path * "mcmcGcov_c.txt", mcmcGcov_c)
    writedlm(analysis_path * "mcmcGcov_total.txt", mcmcGcov_total)

    # save mcmcGcor_c and mcmcGcor_total
    writedlm(analysis_path * "mcmcGcor_c.txt", mcmcGcor_c)
    writedlm(analysis_path * "mcmcGcor_total.txt", mcmcGcor_total)

    # save mcmcAtruecor_c
    writedlm(analysis_path * "mcmcAtruecor_c.txt", mcmcAtruecor_c)

    writedlm(analysis_path * "estAcor.txt", meanAcor)
    writedlm(analysis_path * "estBcor.txt", meanBcor)
    
    if estimate_vara == true
        writedlm(analysis_path * "estBcor_std.txt", sqrt.(abs.(meanBcor2 .- (meanBcor .^ 2))))
        writedlm(analysis_path * "estAcor_std.txt", sqrt.(abs.(meanAcor2 .- (meanAcor .^ 2))))
    end

    # genetic variance components
    for cat in 1:nCategory
        writedlm(analysis_path * "estG" * string(cat) * ".txt", meanG[cat])
        writedlm(analysis_path * "estG_std" * string(cat) * ".txt", sqrt.(abs.(meanG2[cat] .- (meanG[cat] .^ 2))))
        meanGcor[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])])
        meanGcor2[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])] .^ 2)
    end
    writedlm(analysis_path * "estGcor.txt", meanGcor)
    writedlm(analysis_path * "estGcor_std.txt", sqrt.(meanGcor2 .- (meanGcor .^ 2)))

    #  total genetic variance
    writedlm(analysis_path * "estGtotal.txt", meanGtotal, ',')
    writedlm(analysis_path * "estGtotal_std.txt", sqrt.((meanGtotal2 .- (meanGtotal .^ 2))), ',')
    meanGcor_total = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)])
    meanGcor_total2 = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)] .^ 2)
    writedlm(analysis_path * "estGcor_total.txt", meanGcor_total)
    writedlm(analysis_path * "estGcor_total_std.txt", sqrt(meanGcor_total2 - (meanGcor_total^2)))

    if estimate_vare == true
        writedlm(analysis_path * "estR.txt", meanR)
        writedlm(analysis_path * "estR_std.txt", sqrt.((meanR2 .- (meanR .^ 2))))
    end
    if estimate_pi == true
        for cat in 1:nCategory
            SBayesAPP.write_pi_dict(analysis_path * "estPi" * string(cat) * ".txt", mean_pi[cat])
            std_pi = deepcopy(mean_pi[cat])
            for key in SBayesAPP.pi_key_order()
                std_pi[key] = sqrt(mean_pi2[cat][key] - mean_pi[cat][key]^2)
            end
            SBayesAPP.write_pi_dict(analysis_path * "estPi_std" * string(cat) * ".txt", std_pi)
        end
    end

    writedlm(analysis_path * "mcmc_Delta1.rank$my_rank.txt", mcmc_Delta[1])
    writedlm(analysis_path * "mcmc_Delta2.rank$my_rank.txt", mcmc_Delta[2])
    writedlm(analysis_path * "meanAlpha1.rank$my_rank.txt", meanAlpha[1])
    writedlm(analysis_path * "meanAlpha2.rank$my_rank.txt", meanAlpha[2])

    time_end = now()
    time_diff = (time_end - time_start).value / 60000 #milliseconds to min
    println("End time: ", time_end)
    println("Running Time (min): ", time_diff)
run_nonmpi(parse_nonmpi_args(ARGS))
