function parse_metadata(file)
    lines = []
    open(file) do f
        lines = readlines(f)
    end
    parameter_array = split.(lines, "=") |> x->Array{String, 1}.(x)
    parameter_dict = Dict(parameter_array)
    return parameter_dict
end