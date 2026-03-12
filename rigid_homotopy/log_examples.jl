include("rigid_homotopy.jl")
using TaylorSeries
displayBigO(false)
n = 2
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


x1, x2 = set_variables("x1 x2", order=4);

F_taylor = [log(x1+2) log(x2+2)]


eR, vars = make_ring(n)
bind_vars!(vars, n)
HR, (t) = eR["t"]
vars_F = gens(eR)[1:n]
X1 = 0.6931471805599453 + 0.5* x1 - 0.125* x1^2 + 0.041666666666666664* x1^3 - 0.015625* x1^4
X2 = 0.6931471805599453 + 0.5* x2 - 0.125* x2^2 + 0.041666666666666664* x2^3 - 0.015625* x2^4

F = [x1^2 + x2^2 - 1/2 * x1^4-x1^2*x2^2-1/2*x2^4 1.0986122886681098+1.0* x1^2 + 0.3333333333333333 *x1 *x2 - 0.5 *x1^4 - 0.3333333333333333 *x1^3 *x2 - 0.05555555555555555* x1^2 *x2^2]
F = [X1 + X2 X1 - 2*X2]

points = toCCi([[5.08953, 5.08953], [-3.05422, 2-2.54063*im], [7.05422, 2-2.54063*im], [-3.05422,
      2+2.54063*im], [-1.08953, 5.08953], [7.05422, 2+2.54063*im], [5.08953, -1.08953], [-1.08953,
      -1.08953]],CCi)


function rigid_transformations(n)
    transformation_list = []
    for i in 1:n
        M = rand(ComplexF64, n, n)
#        M = Diagonal(map(i -> exp(im*rand(Int)*π/2), 1:n))
        A = matrix(convert_to_box_matrix(M,CCi))
#        A = random_rotation_matrix(n)
#        A = matrix(convert_to_box_matrix([A zeros(ComplexF64,n); zeros(ComplexF64,1,n+1)],CCi))
        push!(transformation_list, matrix_exponential(A, t; order = 4))
    end
    return transformation_list
end
wt = rigid_transformations(n)

iters = 0
pts = []
for i in 1:length(points)
    p, it = rigid_track(F, points[i], .1, wt)
    println(it)
    push!(pts, p)
    iters += it
end
pts, iters

p1, iters1 = rigid_track(F, points[1], .1, wt)
p2, iters2 = rigid_track(F, points[2], .1, wt)



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

