using DelimitedFiles: readdlm
using JLD2

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