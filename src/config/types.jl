module ConfigTypes

export NonMPIConfig, MPIConfig

struct NonMPIConfig
    data_path::String
    analysis_path::String
    nIter::Int
    seed::Int
    nrank::Int
    annot_file::String
    annot_dict::String
    out_freq::Int
    starting_value_dir::String
    secondary_starting_value_dir::String
    st_path::String
    thin::Int
    is_continue::Bool
end

struct MPIConfig
    data_path::String
    analysis_path::String
    nIter::Int
    seed::Int
    nrank::Int
    annot_file::String
    annot_dict::String
    out_freq::Int
    starting_value_dir::String
    secondary_starting_value_dir::String
    st_path::String
    thin::Int
    n1::Int
    n2::Int
    estimate_pi::Bool
    fixed_hyperparameters::Bool
    is_continue::Bool
    chr::String
end

end