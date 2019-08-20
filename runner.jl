using PyPlot # for tests
import Random

include("lobpcg.jl")

N = 500
M = 10
figure()
offset = 1
Random.seed!(offset)
α = 1. # use this to make the problem easier (α < 1) or harder (α > 1)
A = Diagonal(sort(abs.(rand(N)).^α))

β = 4
D = Diagonal(sort(abs.(randn(N)).^β))
B = randn(N,N)
B,_ = qr(B)
B = Array(B)
B = B*D*B'


X,_ = qr(randn(N,M))
X = Array(X)
while norm(X'*B*X - I) > 1e-10
    O = Hermitian(X'*B*X)
    U = cholesky(O).U
    X .= X/U
end

# loose tolerance
tol = 1e-4
ortho_tol = 1e-6 # This uses less orthogonalizations steps than the default
X, resids = LOBPCG(A, X, B, I, tol, 200, ortho_tol=ortho_tol)
figure()
semilogy(resids', "-x")


# # tight tolerance
# tol = 1e-12
# ortho_tol = 2eps()
# X, resids = LOBPCG(A, X, B, I, tol, 200, ortho_tol=ortho_tol)
# figure()
# semilogy(resids', "-x")
