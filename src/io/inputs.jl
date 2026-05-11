using CSV
using DataFrames
using DelimitedFiles: readdlm
using JLD2

function load_annotation_metadata(data_path::AbstractString, annot_file::AbstractString; nCon::Int=0)
    annot = CSV.read(data_path * annot_file, DataFrame)
    annotation_name = String.(names(annot)[2:end])
    n_annotation = length(annotation_name)
    0 <= nCon <= n_annotation || error("nCon must be between 0 and $n_annotation for $(data_path * annot_file)")

    annotation_values = Matrix(annot[!, 2:end])
    n_loci_annot = Int.(vec(sum(annotation_values .!= 0, dims=1)))
    n_cat = n_annotation - nCon
    annotation_type = vcat(fill("continue", nCon), fill("category", n_cat))

    return (
        annotationName=annotation_name,
        nLoci_annot=n_loci_annot,
        nCon=nCon,
        nCat=n_cat,
        annotationType=annotation_type,
    )
end

function load_nonmpi_block_data(data_path::AbstractString, annot_dict::AbstractString)
    transformed_x_dict = JLD2.load(data_path * "TransformedX_dict.jld2")["my_TransformedX_dict"]
    transformed_y_dict = JLD2.load(data_path * "TransformedY_dict.jld2")["my_TransformedY_dict"]
    blkSNPsIndex_dict = JLD2.load(data_path * "blkSNPsIndex_dict.jld2")["my_blkSNPsIndex_dict"]
    blkID = Int.(vec(readdlm(data_path * "blkIDs.txt", ',')))
    nGWAS_dict = JLD2.load(data_path * "nGWAS_dict.jld2")["my_nGWAS_dict"]
    anno_matrix_dict = JLD2.load(data_path * "$annot_dict.jld2")["my_anno_matrix_dict"]
    nblk = length(blkID)
    nsnp = sum(map(length, values(blkSNPsIndex_dict)))

    return (
        transformed_x_dict=transformed_x_dict,
        transformed_y_dict=transformed_y_dict,
        blkSNPsIndex_dict=blkSNPsIndex_dict,
        blkID=blkID,
        nGWAS_dict=nGWAS_dict,
        anno_matrix_dict=anno_matrix_dict,
        nblk=nblk,
        nsnp=nsnp,
    )
end

function load_effect_state(effect_starting_path, delta_starting_path, my_rank, nTraits)
    betaArray = [
        vec(readdlm(effect_starting_path * "last_mcmc_betaArray$(trait).rank$my_rank.txt"))
        for trait in 1:nTraits
    ]
    deltaArray = [
        vec(readdlm(delta_starting_path * "last_sample_delta$(trait)_rank$my_rank.txt"))
        for trait in 1:nTraits
    ]
    alphaArray = [deltaArray[trait] .* betaArray[trait] for trait in 1:nTraits]
    return betaArray, alphaArray, deltaArray
end