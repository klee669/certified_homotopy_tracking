include("rigid_homotopy.jl")
include("../src/certified_monodromy_computation.jl")

n = 3
CCi = AcbField()

function make_ring(n)
    CCi = AcbField()
    varnames = [string("x", i) for i in 1:n]
    append!(varnames, ["λ", "η"])
    return CCi[varnames...]
end
function bind_vars!(vars::Vector, n::Int)
    for i in 1:n
        @eval $(Symbol("x$i")) = vars[$i]
    end
    @eval λ = vars[$n + 1]
    @eval η = vars[$n + 2]
end

eR, vars = make_ring(n)
bind_vars!(vars, n)
HR, (t) = eR["t"]


A_0 = rand(Float64, n, n)
evals0 = eigvals(A_0)
evecs0 = eigvecs(A_0)
A_0 = matrix(convert_to_box_matrix(A_0,CCi))

A_1 = rand(Float64, n, n)
evals1 = eigvals(A_1)
evecs1 = eigvecs(A_1)
points = map(i -> vec(convert_to_box_matrix(hcat([evecs0[:,i]; evals0[i]]), CCi)), 1:n)
A_1 = matrix(convert_to_box_matrix(A_1,CCi))


vars_F = gens(eR)[1:n]
F = (A_0+ λ^2*A_1)*vars_F-vec(Matrix(λ*matrix(vars_F)))
F_homog = (η^2* A_0+ λ^2*A_1)*vars_F-vec(Matrix(η* λ*matrix(vars_F)))

F = Matrix(transpose([F; sum(map(i -> i^2, vars_F))-1]))
F_homog = Matrix(transpose([F_homog; sum(map(i -> i^2, vars_F))-η^2]))


u = map(i -> sample_point(F[i], n+1), 1:n+1)
g = map(i -> gradient_vector(F[i], u[i]), 1:n+1)
frames = map(i -> source_frame_homog(u[i], g[i]), 1:n+1)

zeta = map(i -> CCi(i), rand(ComplexF64, n+1))
target_frames = [];
for i in 1:n+1
    a = map(i -> CCi(i), rand(ComplexF64, n+1))
    b = CCi(rand(ComplexF64))
    
    push!(target_frames, target_frame(zeta,a,b))
end
As = map(i -> log(target_frames[i]*frames[i]'), 1:n+1)
T = sqrt(sum(i -> (1/sqrt(2)*norm(i,2))^2, As))
wt= map(i -> matrix(convert_to_box_matrix(frames[i]*target_frames[i]',CCi))*matrix_exponential(matrix(convert_to_box_matrix(As[i], CCi)), t/T), 1:n+1)
transformations = map(i ->Matrix(i)*vars, wt)
H = map(i -> AbstractAlgebra.evaluate(F_homog[i], transformations[i]), 1:n+1)

point = zeta
r= .1
G = evaluate_matrix(evaluate_matrix(hcat(H), 0), [zeta; CCi(1)])
jacobian_inverse(G, [zeta; CCi(1)])
pseudo_inv(evaluate_matrix(jac(evaluate_matrix(hcat(H), 0)), [zeta; CCi(1)]))












zzz = map(i -> evaluate_matrix(Matrix(wt[i]), T)*[zetas[i]; CCi(1)], 1:n+1)
zz = map(i -> (1/(i[5])*matrix(i)), zzz)


function random_rotation_matrix(n)
    Q, R = qr(randn(n, n))
    Q *= det(Q) < 0 ? -1 : 1  # Make sure it's in SO(n), not just O(n)
    return Matrix(Q)
end
function rigid_transformations(n)
    transformation_list = []
    for i in 1:n
        M = rand(Float64, n, n)
        M = [M zeros(Float64,n); zeros(Float64,1,n+1)]
        A = matrix(convert_to_box_matrix((M- transpose(M))/2,CCi))
#        A = random_rotation_matrix(n)
#        A = matrix(convert_to_box_matrix([A zeros(ComplexF64,n); zeros(ComplexF64,1,n+1)],CCi))
        push!(transformation_list, matrix_exponential(A, t; order = 10))
    end
    return transformation_list
end


wt = rigid_transformations(n)
id_transformation = matrix_exponential(matrix(Matrix(1*I,n+1,n+1)), t; order = 10)
wt = push!(wt, id_transformation)


    R = parent(F[1])
    n = length(F)
    vars_F = gens(R)[1:n]

    w_t = wt
    transformations = map(i -> i*vars_F, w_t)
    transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; CCi(0)]), 1:n)

H = transformed_F
    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
evaluate_matrix(G,points[1])

track(H, points[1], .1; iterations_count = true)
track(H, points[2], .1; iterations_count = true)
track(H, points[3], .1; iterations_count = true)

    R = parent(F[1])
n = length(F)
vars_F = gens(R)[1:n]

w_t = map(i -> evaluate_matrix(Matrix(i), 1), wt)
transformations = map(i -> i*vars_F, w_t)
target_system = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)

HH = vec(Matrix((1-t)*CCi(exp(im*rand(Int)))*matrix(F)+ t*matrix(transpose(target_system))))
lp1, liter1 = track(HH, points[1], .1; iterations_count = true)
lp2, liter2 = track(HH, points[2], .1; iterations_count = true)
lp3, liter3 = track(HH, points[3], .1; iterations_count = true)






# a two-variable example
CCi = AcbField()
eR, (x,y,η) = CCi["x","y","η"]
HR, (t) = eR["t"]

F = [-x^2+y ,  29/16*x^3-2*x*y] + [(1-t)*2, (1-t)*3.8125]
point = [CCi(-1), CCi(-1)]

track(F, point, .1; iterations_count = true)

function rigid_transformations(n)
    transformation_list = []
    for i in 1:n
        M = rand(Float64, n, n)
#        M = [M ; zeros(Float64,1,n+1)]
        A = matrix(convert_to_box_matrix((M- transpose(M))/2,CCi))
        e = exp(-A)
#        A = random_rotation_matrix(n)
#        A = matrix(convert_to_box_matrix([A zeros(ComplexF64,n); zeros(ComplexF64,1,n+1)],CCi))
        push!(transformation_list, e*matrix_exponential(A, (1-t); order = 10))
    end
    return transformation_list
end

n=2 
w_t = rigid_transformations(n)

    R = parent(F[1])
    n = length(F)
    vars_F = gens(eR)[1:n]

    F = [-x^2+y,  29/16*x^3-2*x*y]
    transformations = map(i -> i*vars_F, w_t)
    transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; CCi(0)]), 1:n)
    transformations0 = map(i -> i*point, w_t)
    transformed_F0 = map(i -> AbstractAlgebra.evaluate(F[i], [transformations0[i]; CCi(0)]), 1:n)


H = transformed_F - vec(Matrix((1-t)*matrix(transformed_F0)))
    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
evaluate_matrix(G,point)

track(H, point, .1; iterations_count = true)





# a two-variable example
CCi = AcbField()
eR, (x,y,η) = CCi["x","y","η"]
HR, (t) = eR["t"]

F = [-x^2+y ,  29/16*x^3-2*x*y] + [(1-t)*2, (1-t)*3.8125]
point = [CCi(-1), CCi(-1)]

track(F, point, .1; iterations_count = true)

    n=2
    transformation_list = []
    M1 = [.991434110852759010 -0.13052619222005157; 0.13052619222005157 0.991434110852759010]
    M2 = [0.99622265304372064 -0.08682408883346515; 0.08682408883346515 0.99622265304372064]
    A1 = matrix(convert_to_box_matrix((M1- transpose(M1))/2,CCi))
    A2 = matrix(convert_to_box_matrix((M2- transpose(M2))/2,CCi))
    e1 = exp(-A1)
    e2 = exp(-A2)
    w_t = [e1*matrix_exponential(A1, (1-t); order = 10),
           e2*matrix_exponential(A2, (1-t); order = 10)]

    R = parent(F[1])
    n = length(F)
    vars_F = gens(eR)[1:n]

    F = [-x^2+y,  29/16*x^3-2*x*y]
    transformations = map(i -> i*vars_F, w_t)
    transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; CCi(0)]), 1:n)
    transformations0 = map(i -> i*point, w_t)
    transformed_F0 = map(i -> AbstractAlgebra.evaluate(F[i], [transformations0[i]; CCi(0)]), 1:n)


H = transformed_F - vec(Matrix((1-t)*matrix(transformed_F0)))
    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
evaluate_matrix(G,point)

track(H, point, .1; iterations_count = true)



M1 = [3/5 -4/5; 4/5 3/5]
M2 = [5/13 -12/13;12/13 5/13]
    A1 = matrix(convert_to_box_matrix((M1- transpose(M1))/2,CCi))
    A2 = matrix(convert_to_box_matrix((M2- transpose(M2))/2,CCi))
    e1 = exp(-A1)
    e2 = exp(-A2)
    w_t = [e1*matrix_exponential(A1, (1-t); order = 10),
           e2*matrix_exponential(A2, (1-t); order = 10)]

    R = parent(F[1])
    n = length(F)
    vars_F = gens(eR)[1:n]

    F = [-x^2+y,  29/16*x^3-2*x*y]
    transformations = map(i -> i*vars_F, w_t)
    transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; CCi(0)]), 1:n)
    transformations0 = map(i -> i*point, w_t)
    transformed_F0 = map(i -> AbstractAlgebra.evaluate(F[i], [transformations0[i]; CCi(0)]), 1:n)


H = transformed_F - vec(Matrix((1-t)*matrix(transformed_F0)))
    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
evaluate_matrix(G,point)

track(H, point, .1; iterations_count = true)



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
#    @eval λ = vars[$n + 1]
    @eval η = vars[$n + 1]
end
n=2
eR, vars = make_ring(n)
bind_vars!(vars, n)
HR, (t) = eR["t"]


A_00 = rand(Float64, n, n)
A_00 = matrix(convert_to_box_matrix(A_00,CCi))
A_01 = rand(Float64, n, n)
A_01 = matrix(convert_to_box_matrix(A_01,CCi))
A_02 = rand(Float64, n, n)
A_02 = matrix(convert_to_box_matrix(A_02,CCi))

A_10 = rand(Float64, n, n)
A_10 = matrix(convert_to_box_matrix(A_10,CCi))
A_11 = rand(Float64, n, n)
A_11 = matrix(convert_to_box_matrix(A_11,CCi))
A_12 = rand(Float64, n, n)
A_12 = matrix(convert_to_box_matrix(A_12,CCi))



vars_F = gens(eR)[1:n]
F = [det(A_00*x+A_01*y+A_02); det(A_10*x+A_11*y+A_12)]
#F = Matrix(transpose([F; sum(map(i -> i^2, vars_F))-1]))



function random_rotation_matrix(n)
    Q, R = qr(randn(n, n))
    Q *= det(Q) < 0 ? -1 : 1  # Make sure it's in SO(n), not just O(n)
    return Matrix(Q)
end
function rigid_transformations(n)
    transformation_list = []
    for i in 1:n
        M = rand(Float64, n, n)
#        M = [M zeros(Float64,n); zeros(Float64,1,n+1)]
        A = matrix(convert_to_box_matrix((M- transpose(M))/2,CCi))
#        A = random_rotation_matrix(n)
#        A = matrix(convert_to_box_matrix([A zeros(ComplexF64,n); zeros(ComplexF64,1,n+1)],CCi))
        push!(transformation_list, matrix_exponential(A, t; order = 10))
    end
    return transformation_list
end


wt = rigid_transformations(n)


    R = parent(F[1])
    n = length(F)
    vars_F = gens(R)[1:n]


    transformations = map(i -> i*vars_F, w_t)
    transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; CCi(0)]), 1:n)

H = transformed_F
    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
evaluate_matrix(G,points[1])

track(H, points[1], .1; iterations_count = true)
track(H, points[2], .1; iterations_count = true)
track(H, points[3], .1; iterations_count = true)