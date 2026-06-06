using DelimitedFiles,DataFrames,CSV,Statistics
  
h2_all = [(0.5, 0.2)]

res_folder_path = "/common/zhao/jyqqu/MTSBayesCC/analysis/sim_chr1_output_v2/MTS/MTSBayesCC_SNPorder_ShuffleTraits/"
data_folder_path = "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/"
save_path = "/common/zhao/jyqqu/MTSBayesCC/analysis/sim_chr1_output_v2/MTS/"

nind_col     = []
plepct_col   = []
seed_col     = []

gcov_c1_simu_col = []
gcov_c2_simu_col = []

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

                    for c in 1:2
                        bv_var_file_path = data_path * "bvc$(c)_cov_var.txt"
                        saved_bv_var = readdlm(bv_var_file_path)
                        gvar_simu[c] = saved_bv_var
                    end

		    append!(nind_col,samSize)
                    append!(plepct_col,plePer)
                    append!(seed_col,seed)

		   append!(gcov_c1_simu_col, gvar_simu[1][1,2])
                   append!(gcov_c2_simu_col, gvar_simu[2][1,2])
		end		
	   end
	end
	df=DataFrame(nind=nind_col, plepct=plepct_col, seed=seed_col,
                gcov_c1_simu=gcov_c1_simu_col,gcov_c2_simu=gcov_c2_simu_col)
	CSV.write(save_path*"anno_stratified_true_gcov_h2.h2_trait1.$h21.h2_trait2.$h22.txt",df,delim='\t')
end


