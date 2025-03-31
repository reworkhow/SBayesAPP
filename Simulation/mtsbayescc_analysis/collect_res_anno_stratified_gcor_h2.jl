using DelimitedFiles,DataFrames,CSV,Statistics


data_folder_path = "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/"
res_folder_path = "/common/zhao/jyqqu/MTSBayesCC/analysis/sim_chr1_output_v2/MTS/MTSBayesCC_SNPorder_ShuffleTraits/" #MT-S-BayesC-C (EIEO)
save_path = "$res_folder_path/res_summary/"
#mkpath(save_path)

# aiao_folder_path = "/common/zhao/tianjing/simu_ch1_res/mt_s_bayesc_c_aiao/" #MT-S-BayesC-C (AIAO)
# mt_s_bayesc_folder_path = "/common/zhao/tianjing/simu_ch1_res/mt_s_bayesc/" #MT-S-BayesC (no annot)
# rrblup_folder_path = "/common/zhao/tianjing/simu_ch1_res/mt_s_rrblup/" #MT-S-RRBLUP (no annot)


h2_all = [(0.5, 0.2)]

nind_col     = []
#annosize_col = []
plepct_col   = []
seed_col     = []

gcor_c1_simu_col = []
gcor_c2_simu_col = []

marker_cor_c1_simu_col = []
marker_cor_c2_simu_col = []

h21_c1_simu_col = []
h22_c1_simu_col = []
h21_c2_simu_col = []
h22_c2_simu_col = []

gcor_c1_est_col = []
gcor_c2_est_col = []

marker_cor_c1_est_col = []
marker_cor_c2_est_col = []

h21_c1_est_col = []
h22_c1_est_col = []
h21_c2_est_col = []
h22_c2_est_col = []

gcor_c1_gnova_col=[]
gcor_c2_gnova_col=[]

h21_c1_gnova_col=[]
h22_c1_gnova_col=[]
h21_c2_gnova_col=[]
h22_c2_gnova_col=[]

# gcor_c1_aiao_col=[]
# gcor_c2_aiao_col=[]

for (h21, h22) in h2_all
    for samSize in [300000]
        #for annoSize in [0.5]
            for plePer in [0.1,0.5,0.8,1]
                if plePer==0 || plePer==1
                    plePer=Int64(plePer)
                end
                for seed in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
                    #get varg from data
                    folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(plePer)_sampleSize$(samSize)_seed$seed/"
                    println(folder_name)
                    data_path = data_folder_path*folder_name
                    
                    gvar_simu = [zeros(2, 2) for c in 1:2]
		    marker_var_simu = [zeros(2, 2) for c in 1:2]

                    # saved_sdbv = readdlm(data_path * "sdbv.txt")[:, 1]
                    # saved_sdy = readdlm(data_path * "sdy.txt")[:, 1]

                    for c in 1:2
                        bv_var_file_path = data_path * "bvc$(c)_cov_var.txt"
			marker_var_file_path = data_path * "QTLcovmatc$(c).txt"
                        saved_bv_var = readdlm(bv_var_file_path)
			saved_marker_var = readdlm(marker_var_file_path)
                        gvar_simu[c] = saved_bv_var  
			marker_var_simu[c] = saved_marker_var
                    end

                    #estimated genetic correlation from mt-s-bayesc-c eieo
                    gcor_res_path = res_folder_path * folder_name*"/"*"estGcor.txt"
		    marker_cor_res_path = res_folder_path * folder_name*"/"*"estBcor.txt"

                    h2_c1_res_path = res_folder_path * folder_name * "/" * "estG1.txt"
                    h2_c2_res_path = res_folder_path * folder_name * "/" * "estG2.txt"

                    if isfile(gcor_res_path)
                        gcor_c1c2_est = vec(readdlm(gcor_res_path))
			marker_cor_c1c2_est = vec(readdlm(marker_cor_res_path))
                        gcor_c1_est = gcor_c1c2_est[1]
                        gcor_c2_est = gcor_c1c2_est[2]
			marker_cor_c1_est = marker_cor_c1c2_est[1]
			marker_cor_c2_est = marker_cor_c1c2_est[2]
                        h21_c1_est = readdlm(h2_c1_res_path)[1, 1]
                        h22_c1_est = readdlm(h2_c1_res_path)[2, 2]
                        h21_c2_est = readdlm(h2_c2_res_path)[1, 1]
                        h22_c2_est = readdlm(h2_c2_res_path)[2, 2]
                    else
                        @warn "Estimate File not found: $folder_name"
                        gcor_c1_est = gcor_c2_est =  999
			marker_c1_est = marker_c2_est = 999
                        h21_c1_est = h22_c1_est = 999
                        h21_c2_est = h22_c2_est = 999
                    end
                    

                    #estimated genetic correlation from mt-s-bayesc-c aiao
                    # res_path2   = aiao_folder_path*folder_name*"/"*"estGcor.txt"
                    # gcor_c1c2_aiao = vec(readdlm(res_path2))
                    # gcor_c1_aiao = gcor_c1c2_aiao[1]
                    # gcor_c2_aiao = gcor_c1c2_aiao[2]
                    
                    #estimated genetic correlation from gnova
                    gnova_res_path = data_path * "/" * "gnova.txt"
                    gnova_res = CSV.read(gnova_res_path, DataFrame)
                    gcor_c1_gnova = gnova_res[1, :corr_corrected]
                    gcor_c2_gnova = gnova_res[2, :corr_corrected]
                    h21_c1_gnova = gnova_res[1, :h2_1]
                    h21_c2_gnova = gnova_res[2, :h2_1]
                    h22_c1_gnova = gnova_res[1, :h2_2]
                    h22_c2_gnova = gnova_res[2, :h2_2]

                    #gnova may have NA
                    if gcor_c1_gnova =="NA"
                        gcor_c1_gnova = 999
                    end
                    if gcor_c2_gnova =="NA"
                        gcor_c2_gnova = 999
                    end
                    if typeof(gcor_c1_gnova) <: AbstractString
                        gcor_c1_gnova = parse(Float64, gcor_c1_gnova)
                    end
                    if typeof(gcor_c2_gnova) <: AbstractString
                        gcor_c2_gnova = parse(Float64, gcor_c2_gnova)
                    end
                    if typeof(h21_c1_gnova) <: AbstractString
                        h21_c1_gnova = parse(Float64, h21_c1_gnova)
                    end
                    if typeof(h21_c2_gnova) <: AbstractString
                        h21_c2_gnova = parse(Float64, h21_c2_gnova)
                    end
                    if typeof(h22_c1_gnova) <: AbstractString
                        h22_c1_gnova = parse(Float64, h22_c1_gnova)
                    end
                    if typeof(h22_c2_gnova) <: AbstractString
                        h22_c2_gnova = parse(Float64, h22_c2_gnova)
                    end

                    append!(nind_col,samSize)
                    #append!(annosize_col,annoSize)
                    append!(plepct_col,plePer)
                    append!(seed_col,seed)

                    append!(gcor_c1_simu_col, gvar_simu[1][1,2]/sqrt(gvar_simu[1][1,1]*gvar_simu[1][2,2]))
                    append!(gcor_c2_simu_col, gvar_simu[2][1,2]/sqrt(gvar_simu[2][1,1]*gvar_simu[2][2,2]))
                    append!(marker_cor_c1_simu_col, marker_var_simu[1][1,2]/sqrt(marker_var_simu[1][1,1] * marker_var_simu[1][2,2]))
		    append!(marker_cor_c2_simu_col, marker_var_simu[2][1,2]/sqrt(marker_var_simu[2][1,1] * marker_var_simu[2][2,2]))

		    append!(h21_c1_simu_col, gvar_simu[1][1,1])
                    append!(h22_c1_simu_col, gvar_simu[1][2,2])
                    append!(h21_c2_simu_col, gvar_simu[2][1,1])
                    append!(h22_c2_simu_col, gvar_simu[2][2,2])

                    append!(gcor_c1_est_col,  gcor_c1_est)
                    append!(gcor_c2_est_col,  gcor_c2_est)
		    append!(marker_cor_c1_est_col,  marker_cor_c1_est)
		    append!(marker_cor_c2_est_col,  marker_cor_c2_est)
                    append!(h21_c1_est_col,  h21_c1_est)
                    append!(h22_c1_est_col,  h22_c1_est)
                    append!(h21_c2_est_col,  h21_c2_est)
                    append!(h22_c2_est_col,  h22_c2_est)

                    append!(gcor_c1_gnova_col,  gcor_c1_gnova)
                    append!(gcor_c2_gnova_col,  gcor_c2_gnova)
                    append!(h21_c1_gnova_col,  h21_c1_gnova)
                    append!(h22_c1_gnova_col,  h22_c1_gnova)
                    append!(h21_c2_gnova_col,  h21_c2_gnova)
                    append!(h22_c2_gnova_col,  h22_c2_gnova)

                    # append!(gcor_c1_aiao_col,  gcor_c1_aiao)
                    # append!(gcor_c2_aiao_col,  gcor_c2_aiao)

                end
            end
        #end
    end
    df=DataFrame(nind=nind_col, plepct=plepct_col, seed=seed_col,
                gcor_c1_simu=gcor_c1_simu_col,gcor_c2_simu=gcor_c2_simu_col,
                h21_c1_simu=h21_c1_simu_col,h22_c1_simu=h22_c1_simu_col,
                h21_c2_simu=h21_c2_simu_col,h22_c2_simu=h22_c2_simu_col,
                gcor_c1_est=gcor_c1_est_col, gcor_c2_est=gcor_c2_est_col,
                h21_c1_est=h21_c1_est_col, h22_c1_est=h22_c1_est_col,
                h21_c2_est=h21_c2_est_col, h22_c2_est=h22_c2_est_col,
                gcor_c1_gnova=gcor_c1_gnova_col, gcor_c2_gnova=gcor_c2_gnova_col,
                h21_c1_gnova=h21_c1_gnova_col, h22_c1_gnova=h22_c1_gnova_col,
                h21_c2_gnova=h21_c2_gnova_col, h22_c2_gnova=h22_c2_gnova_col,
		marker_cor_c1_est=marker_cor_c1_est_col, marker_cor_c2_est=marker_cor_c2_est_col,
		marker_cor_c1_simu=marker_cor_c1_simu_col,marker_cor_c2_simu=marker_cor_c2_simu_col)
                # gcor_c1_aiao=gcor_c1_aiao_col, gcor_c2_aiao=gcor_c2_aiao_col
                # annosize=annosize_col,

    CSV.write(save_path*"anno_stratified_gcor_h2.h2_trait1.$h21.h2_trait2.$h22.txt",df,delim='\t')
end









