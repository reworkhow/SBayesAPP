using DelimitedFiles: readdlm

function build_start_pi(st_path::AbstractString; estimate_pi::Bool=true)
    if estimate_pi
        pi11 = 0.00001
        pi10 = 1.0 - round(readdlm(st_path * "Trait1/mean_pi.txt")[1, 1], digits=4)
        pi01 = 1.0 - round(readdlm(st_path * "Trait2/mean_pi.txt")[1, 1], digits=4)
        pi00 = 1.0 - pi11 - pi10 - pi01
        start_pi = Dict([1.0; 1.0] => pi11, [1.0; 0.0] => pi10, [0.0; 1.0] => pi01, [0.0; 0.0] => pi00)
    else
        pi00 = 1e-04
        start_pi = Dict([1.0; 1.0] => 0.9997, [1.0; 0.0] => 1e-04, [0.0; 1.0] => 1e-04, [0.0; 0.0] => 1e-04)
    end

    return (startPi=start_pi, Pi00=pi00)
end

function build_gprior_vec(st_path::AbstractString, n_loci_annot, pi00)
    st_h21 = round(readdlm(st_path * "Trait1/mean_varg_total.txt")[1, 1], digits=3)
    st_h22 = round(readdlm(st_path * "Trait2/mean_varg_total.txt")[1, 1], digits=3)

    gprior_vec = [zeros(2, 2) for _ in 1:length(n_loci_annot)]
    loci_total = sum(n_loci_annot)
    for category in eachindex(n_loci_annot)
        gprior_vec[category] = [st_h21 0.0; 0.0 st_h22] * (n_loci_annot[category] / loci_total)
        gprior_vec[category] = gprior_vec[category] / (n_loci_annot[category] * (1 - pi00))
    end

    return gprior_vec
end