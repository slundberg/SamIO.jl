module SamIO

abstract SequenceReader

include("ReferenceContigs.jl")
include("BamReader.jl")
include("BedReader.jl")
include("BinningMap.jl")

end # module
