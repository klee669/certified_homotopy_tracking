include("rigid_homotopy.jl")

n = 3
CCi = AcbField()

function make_ring(n)
    CCi = AcbField()
    varnames = [string("x", i) for i in 1:n]
    append!(varnames, ["η"])
    return CCi[varnames...]
end
function bind_vars!(vars::Vector, n::Int)
    for i in 1:n
        @eval $(Symbol("x$i")) = vars[$i]
    end
    @eval η = vars[$n + 1]
end

eR, vars = make_ring(n)
bind_vars!(vars, n)
HR, (t) = eR["t"]


A_0 = rand(ComplexF64, n, n)
evals = eigvals(A_0)
evecs = eigvecs(A_0)
points = map(i -> vec(convert_to_box_matrix(hcat(evecs[:,i]), CCi)), 1:n)
A_0 = matrix(convert_to_box_matrix(A_0,CCi))

A_0 = rand(ComplexF64, n, n)
A_1 = matrix(convert_to_box_matrix(A_0,CCi))

A_0 = rand(ComplexF64, n, n)
A_2 = matrix(convert_to_box_matrix(A_0,CCi))


vars_F = gens(eR)[1:n]

F1 = transpose(matrix(vars_F))*matrix(A_0*vars_F)
F2 = transpose(matrix(vars_F))*matrix(A_1*vars_F)
F3 = transpose(matrix(vars_F))*matrix(A_2*vars_F)
F = Matrix([F1 F2 F3])



function random_rotation_matrix(n)
    Q, R = qr(randn(n, n))
    Q *= det(Q) < 0 ? -1 : 1  # Make sure it's in SO(n), not just O(n)
    return Matrix(Q)
end
function rigid_transformations(n)
    transformation_list = []
    for i in 1:n
        M = rand(ComplexF64, n, n)
#        M = Diagonal(map(i -> exp(im*rand(Int)*π/2), 1:n))
#        M = [M zeros(ComplexF64,n); zeros(ComplexF64,1,n)]
        A = matrix(convert_to_box_matrix((M- M')/2,CCi))
#        A = random_rotation_matrix(n)
#        A = matrix(convert_to_box_matrix([A zeros(ComplexF64,n); zeros(ComplexF64,1,n+1)],CCi))
        push!(transformation_list, matrix_exponential(A, t; order = 10))
    end
    return transformation_list
end


wt = rigid_transformations(n)
id_transformation = matrix_exponential(matrix(Matrix(1*I,n,n)), t; order = 10)
wt = push!(wt, id_transformation)


p1, iters1 = rigid_track(F, points[1], .1, wt)
p2, iters2 = rigid_track(F, points[2], .1, wt)
p3, iters3 = rigid_track(F, points[3], .1, wt)
iters1 + iters2 + iters3

# comparison with the certified linear homotopy
# These are needed to construct the linear homotopy
    R = parent(F[1])
n = length(F)
vars_F = gens(R)[1:n]

w_t = map(i -> evaluate_matrix(Matrix(i), 1), wt)
transformations = map(i -> i*vars_F, w_t)
target_system = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)

H = Matrix((1-t)*matrix(transpose(F))+ t*matrix(target_system))
HH = Matrix((1-t)*matrix(F)+ t*matrix(transpose(target_system)))
lp1, liter1 = tracking_without_predictor(HH, points[1], .1; iterations_count = true)
lp2, liter2 = tracking_without_predictor(HH, points[2], .1; iterations_count = true)
lp3, liter3 = tracking_without_predictor(HH, points[3], .1; iterations_count = true)
liter1 + liter2 + liter3

