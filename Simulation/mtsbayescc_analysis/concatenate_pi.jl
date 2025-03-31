using CSV, DataFrames
res_folder = ARGS[1]

# Function to read estPi files and return a dictionary
function read_estPi(file)
    pi_dict = Dict("[0.0, 0.0]" => NaN, "[1.0, 1.0]" => NaN, "[1.0, 0.0]" => NaN, "[0.0, 1.0]" => NaN)
    if isfile(file)
        for line in eachline(file)
            key, value = split(line, "\t", limit=2)
            if haskey(pi_dict, key)
                pi_dict[key] = parse(Float64, value)
            end
        end
    end
    return pi_dict
end

# Initialize a DataFrame to store the results
results = DataFrame(Folder=String[], Pi00_c1=Float64[], Pi11_c1=Float64[],
    Pi10_c1=Float64[], Pi01_c1=Float64[], Pi00_c2=Float64[],
    Pi11_c2=Float64[], Pi10_c2=Float64[], Pi01_c2=Float64[])

# Iterate over subfolders
for folder in readdir(res_folder, join=true)
    if isdir(folder)
        folder_name = basename(folder)

        # Read values from estPi1.txt and estPi2.txt
        pi1 = read_estPi(joinpath(folder, "estPi1.txt"))
        pi2 = read_estPi(joinpath(folder, "estPi2.txt"))

        # Add a row to the DataFrame
        push!(results, (folder_name, pi1["[0.0, 0.0]"], pi1["[1.0, 1.0]"], pi1["[1.0, 0.0]"], pi1["[0.0, 1.0]"],
            pi2["[0.0, 0.0]"], pi2["[1.0, 1.0]"], pi2["[1.0, 0.0]"], pi2["[0.0, 1.0]"]))
    end
end

# Save the DataFrame to a CSV file
CSV.write("$(res_folder)/pi_values.csv", results)
println("Extraction complete! Results saved to pi_values.csv.")
