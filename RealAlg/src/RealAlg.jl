__precompile__()

module RealAlg

using SparseArrays

import Base: +,-,*,==,>,>=,<,<=,broadcast,sin,cos,tan,cot,abs,exp,log,log10
export Cell C_e

include("RealAlgTypes.jl")
include("Functions/C_e.jl")

end # module