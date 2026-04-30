using CSV
using DataFrames

function load_annotation_metadata(data_path::AbstractString, annot_file::AbstractString; nCon::Int=0)
    annot = CSV.read(data_path * annot_file, DataFrame)
    annotation_name = names(annot)[2:end]
    n_loci_annot = sum.(eachcol(annot[!, 2:end]))
    n_cat = length(annotation_name)
    annotation_type = repeat(["category"], n_cat)

    return (
        annotationName=annotation_name,
        nLoci_annot=n_loci_annot,
        nCon=nCon,
        nCat=n_cat,
        annotationType=annotation_type,
    )
end