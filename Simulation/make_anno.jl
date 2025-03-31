using DataFrames
using DelimitedFiles, CSV


using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed"
        help = ""
        arg_type = Int
        default = 123

        "--pleio_percent"
        help = ""
        arg_type = Float64
        default = 0.1

        "--sample_size"
        help = ""
        arg_type = Int
        default = 10

        # "--annotation_size"
        # help = ""
        # arg_type = Float64
        # default = 0.1

        "--h21"
        help = ""
        arg_type = Float64
        default = 0.01

        "--h22"
        help = ""
        arg_type = Float64
        default = 0.01
    end

    return parse_args(s)
end

# Use the parsed arguments
args = parse_commandline()

# Assuming `args` is the dictionary returned from the parse_args() function
samSize = args["sample_size"]
# annoSize = args["annotation_size"]
plePer = args["pleio_percent"]
if plePer == 0 || plePer == 1
    plePer = Int(plePer) #0, instead of 0.0
end
seed = args["seed"]
h21 = args["h21"] # h2 for trait 1
h22 = args["h22"] # h2 for trait 2

# Optionally, you can print these values to check
println("Sample size: ", samSize)
# println("Annotation size: ", annoSize)
println("Pleiotropy percent: ", plePer)
println("Seed: ", seed)
println("Heritability for trait 1 (h21): ", h21)
println("Heritability for trait 2 (h22): ", h22)
println("----check point a")


basePath = "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.$(h21).h2_trait2.$(h22)_pleioPercent$(plePer)_sampleSize$(samSize)_seed$(seed)/"
snp_list_Path = "/common/zhao/jyqqu/MTSBayesCC/data/bfiles/"
@show basePath
cd(basePath)

snp_c1 = vec(readdlm("../SNPc1.txt"))
snp_c2 = vec(readdlm("../SNPc2.txt"))

snp_list = vec(readdlm(snp_list_Path * "snp_list.1000G.EUR.QC.1.filtered.txt"))
n_snp = length(snp_list)

if length(snp_c1) + length(snp_c2) != n_snp
    error("length not match")
end

c1 = [snp in snp_c1 ? 1 : 0 for snp in snp_list]
c2 = ones(Int, n_snp) - c1

anno = DataFrame(c1=c1, c2=c2)
CSV.write("anno.txt", anno, delim='\t')

