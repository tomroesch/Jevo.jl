
"""
    parse_metadata(file)

Transform a metadata file into a dictionary.

File has to have one parameter per line. Returns a dictionary where parameter names
are stored as keys as strings and the parameter values are stored as values as strings.


# Examples
```julia-repl
shell> cat "test_Metadata.txt"
parameter_1=1
parameter_2=[1, 1]

julia> Jevo.parse_metadata("test_Metadata.txt")
Dict{String,String} with 2 entries:
  "parameter_1"         => "1"
  "parameter_2"         => "[1, 1]"
```
"""
function parse_metadata(file)
    lines = []
    open(file) do f
        lines = readlines(f)
    end
    parameter_array = split.(lines, "=") |> x->Array{String, 1}.(x)
    parameter_dict = Dict(parameter_array)
    return parameter_dict
end