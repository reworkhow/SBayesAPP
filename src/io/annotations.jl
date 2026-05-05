using CSV
using DataFrames

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