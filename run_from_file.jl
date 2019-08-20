using PyPlot
using JLD
using LinearAlgebra
include("lobpcg.jl")

struct DataPrec
    invdata
end
import Base: \
LinearAlgebra.ldiv!(Y, block::DataPrec, B) = mul!(Y, block.invdata, B)
\(block::DataPrec, B) = ldiv!(similar(B), block, B)

function run_from_file(file)
    data = JLD.load(file)
    ortho_tol = max(5eps(real(eltype(data["hamk"]))), data["tol"] / 1000)
    prec = DataPrec(data["inv_preck"])
    X, resids = LOBPCG(data["hamk"], data["guessk"], I, prec, data["tol"], 100;
                       ortho_tol=ortho_tol)

    figure()
    semilogy(resids', "-x")
end
