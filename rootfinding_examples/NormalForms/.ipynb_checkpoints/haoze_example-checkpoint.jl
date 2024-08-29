using MultivariatePolynomials
using DynamicPolynomials
using MAT
include("denseTNF.jl")

#-------------------------------------------------------------------------------
# Intersecting generic curves in ℂ²
n = 3;
d = [10,10,10]; # Intersecting a degree 7 and a degree 13 curve
@polyvar x[1:n]

# Construct a generic system
f = fill(1.0+sum(x),n)
for i = 1:n
    mons = monomials((1+sum(x))^d[i])
    f[i] = randn(length(mons))'*mons;
end

@time sol_QR,M = solveDense(f,x)

file = matopen("commuting_matrics_from_polynomials.mat", "w")
write(file, "M", M)
close(file)