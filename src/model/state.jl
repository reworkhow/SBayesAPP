using DelimitedFiles: readdlm
using LinearAlgebra: inv

function initialize_nonmpi_parameter_state(
    my_rank,
    nCategory,
    nTraits,
    Gprior_vec,
    startPi;
    n1,
    n2,
    estimate_vare=true,
    estimate_vara=true,
    estimate_Gscale=true,
    is_continue=false,
    starting_value_dir="",
    gscale_value_dir="",
)
    A_vec = [zeros(2, 2) for _ in 1:nCategory]
    if is_continue
        A_vec_starting_path = starting_value_dir * "beta_effect_var_matrices_last_sample/"
        for category in 1:nCategory
            A_vec[category] = readdlm(A_vec_starting_path * "beta_effect_matrix_$category.txt", ',')
        end
    else
        A_vec = Gprior_vec
    end
    Ainv_vec = [inv(A_vec[category]) for category in 1:nCategory]

    if is_continue
        Pi_starting_path = starting_value_dir * "pi_last_sample/"
        Pi = [Dict{NTuple{2,Float64},Float64}() for _ in 1:nCategory]
        for category in 1:nCategory
            Pi[category] = read_to_dict(Pi_starting_path * "pi_$category.txt")
        end
    else
        Pi = [deepcopy(startPi) for _ in 1:nCategory]
    end

    Rprior = [1.0 / n1 0.0; 0.0 1.0 / n2]
    df_R = nothing
    scale_R = nothing
    if estimate_vare
        df_R = 4 + nTraits
        scale_R = Rprior * (df_R - nTraits - 1)
    end

    df_G = nothing
    scale_G_vec = nothing
    if my_rank == 0 && estimate_vara
        df_G = 4 + nTraits
        if is_continue && !isempty(gscale_value_dir)
            scale_G_vec = [zeros(2, 2) for _ in 1:nCategory]
            for category in 1:nCategory
                scale_G_vec[category] = readdlm(gscale_value_dir * "scale_G$category.txt")
            end
            estimate_Gscale = false
            println("scale_G_vec is fixed as saved in $gscale_value_dir")
        else
            scale_G_vec = [Gprior_vec[category] * (df_G - nTraits - 1) for category in 1:nCategory]
            println("starting value of scale_G_vec is computed by ST h2.")
        end
    end

    return (
        A_vec=A_vec,
        Ainv_vec=Ainv_vec,
        Pi=Pi,
        Rprior=Rprior,
        df_R=df_R,
        scale_R=scale_R,
        df_G=df_G,
        scale_G_vec=scale_G_vec,
        estimate_Gscale=estimate_Gscale,
    )
end