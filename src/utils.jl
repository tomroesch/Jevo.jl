using InteractiveUtils

export showtypetree
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


"""
    function showtypetree(T, level=0)

Show the all subtypes in a tree. Taken from  https://en.wikibooks.org/wiki/Introducing_Julia/Types .

# Example
```julia-repl
julia> showtypetree(Real)
Real
	AbstractFloat
		BigFloat
		Float16
		Float32
		Float64
	AbstractIrrational
		Irrational
	Integer
		Bool
		Signed
			BigInt
			Int128
			Int16
			Int32
			Int64
			Int8
		Unsigned
			UInt128
			UInt16
			UInt32
			UInt64
			UInt8
	Rational
	StatsBase.TestStat
```
"""
function showtypetree(T, level=0)
    println("\t" ^ level, T)
    for t in subtypes(T)
        showtypetree(t, level+1)
    end
end
