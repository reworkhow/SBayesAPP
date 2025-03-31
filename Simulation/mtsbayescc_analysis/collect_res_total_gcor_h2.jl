using DelimitedFiles,DataFrames,CSV,Statistics

h2_all = [(0.5, 0.2)]

res_folder_path = "/common/zhao/jyqqu/MTSBayesCC/analysis/sim_chr1_output_v2/MTS/MTSBayesCC_SNPorder_ShuffleTraits_estGscale_STinput2/" #MT-S-BayesC-C (EIEO)
data_folder_path = "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/"
save_path = "$res_folder_path/res_summary/"
mkpath(save_path)

# rrblup_folder_path = "/common/zhao/tianjing/simu_ch1_res/mt_s_rrblup/" #MT-S-RRBLUP (no annot)
# aiao_folder_path = "/common/zhao/tianjing/simu_ch1_res/mt_s_bayesc_c_aiao/" #MT-S-BayesC-C (AIAO)
# mt_s_bayesc_folder_path = "/common/zhao/tianjing/simu_ch1_res/mt_s_bayesc/" #MT-S-BayesC (no annot)


nind_col     = []
#annosize_col = []
plepct_col   = []
seed_col     = []

gcor_total_simu_col = []
gcor_total_est_col = []
gcor_ldsc_col = []

# gcor_rrblup_col = []
# gcor_aiao_col = []
# gcor_mt_s_bayesc_col=[]

h21_total_simu_col = []
h21_total_est_col = []
h21_ldsc_col = []

h22_total_simu_col = []
h22_total_est_col = []
h22_ldsc_col = []

ldsc_df = CSV.read(data_folder_path * "summary_ldsc_results.csv", DataFrame)

for (h21, h22) in h2_all
    for samSize in [300000]
        #for annoSize in [0.5]
            for plePer in [0.1,0.5,0.8,1]
                if plePer==0 || plePer==1
                    plePer=Int64(plePer)
                end
                for seed in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
                    ### Total genetic variance
                    #get varg from data 
                    folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(plePer)_sampleSize$(samSize)_seed$seed"
                    println(folder_name)
                    data_path = data_folder_path*folder_name
                    #read bv_trait1 and bv_trait2 to obtain a good starting value for genetic_variance, to make G_scale good.
                    bv_var_file_path = data_path *"/"* "bv_cov_var.txt"
                    if isfile(bv_var_file_path)
                        bv_var = readdlm(bv_var_file_path)
                    else
                        @warn "$bv_var_file_path not found"
                    end
                    gcor_simu = bv_var[1, 2] / sqrt(bv_var[1, 1] * bv_var[2, 2])
                    h21_simu = bv_var[1, 1]
                    h22_simu = bv_var[2, 2]

                    # estimated genetic correlation from mtsbayescc
                    gcor_res_path = res_folder_path * folder_name*"/"*"estGcor_total.txt"
                    h2_res_path = res_folder_path * folder_name * "/" * "estGtotal.txt"
                    if isfile(gcor_res_path)
                        gcor_est = readdlm(gcor_res_path)[1]
                        h21_est = readdlm(h2_res_path, ',')[1,1]
                        h22_est = readdlm(h2_res_path, ',')[2,2]
                    else
                        @warn "Estimate File not found: $folder_name"
                        gcor_est = 999
                        h21_est = 999
                        h22_est = 999
                    end

                    #ldsc
	  	    ldsc_df_i = ldsc_df[ldsc_df[!,:Subfolder] .== folder_name,:]
		    ldsc_gcor = ldsc_df_i[1,:genetic_corr]
		    ldsc_h21 = ldsc_df_i[1,:h2_trait1]
		    ldsc_h22 = ldsc_df_i[1,:h2_trait2]

                   
                    # #estimated genetic correlation from mt-s-rrblup
                    # res_path2   = rrblup_folder_path*folder_name*"/"*"estGcor_total.txt"
                    # gcor_rrblup = readdlm(res_path2)[1]

                    # #estimated genetic correlation from mt_s_bayesc_c_aiao
                    # res_path3   = aiao_folder_path*folder_name*"/"*"estGcor_total.txt"
                    # gcor_aiao = readdlm(res_path3)[1]

                    # #estimated genetic correlation from gcor_mt_s_bayesc
                    # res_path4   = mt_s_bayesc_folder_path*folder_name*"/"*"estGcor_total.txt"
                    # gcor_mt_s_bayesc = readdlm(res_path4)[1]

                    
                    append!(nind_col,samSize)
                    #append!(annosize_col,annoSize)
                    append!(plepct_col,plePer)
                    append!(seed_col,seed)

                    append!(gcor_total_simu_col, gcor_simu)
                    append!(h21_total_simu_col, h21_simu)
                    append!(h22_total_simu_col, h22_simu)
                    
                    append!(gcor_total_est_col,  gcor_est)
                    append!(h21_total_est_col, h21_est)
                    append!(h22_total_est_col, h22_est)

		    append!(gcor_ldsc_col, ldsc_gcor)
		    append!(h21_ldsc_col, ldsc_h21)
		    append!(h22_ldsc_col, ldsc_h22)

                    # append!(gcor_rrblup_col, gcor_rrblup)
                    # append!(gcor_aiao_col, gcor_aiao)
                    # append!(gcor_mt_s_bayesc_col, gcor_mt_s_bayesc)
                end
            end
        #end
    end
    df=DataFrame(nind=nind_col, plepct=plepct_col, seed=seed_col,
             gcor_total_simu=gcor_total_simu_col, h21_total_simu=h21_total_simu_col, h22_total_simu=h22_total_simu_col,
             gcor_total_est=gcor_total_est_col, h21_total_est=h21_total_est_col, h22_total_est=h22_total_est_col,
	     gcor_ldsc=gcor_ldsc_col, h21_ldsc = h21_ldsc_col, h22_ldsc = h22_ldsc_col)
             # annosize=annosize_col, 
             # , gcor_rrblup=gcor_rrblup_col, gcor_aiao=gcor_aiao_col,
             # gcor_mtsbayesc=gcor_mt_s_bayesc_col)

    CSV.write(save_path*"total_gcor_h2.h2_trait1.$h21.h2_trait2.$h22.txt",df,delim='\t')
end











