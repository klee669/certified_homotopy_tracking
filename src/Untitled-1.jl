include("certified_monodromy_computation.jl")


#Katsura 3
CCi = AcbField()
R, (a, b, c, η) = CCi["a","b","c","η"]
HR, (t) = R["t"]

f1 = a+2*b+2*c-1
f2 = a^2+2*b^2+2*c^2-a
f3 =  2*a*b+2*b*c-b

g1 = a-1
g2 = b^2-1
g3 = c^2-1
G = [g1 g2 g3]
F = [f1 f2 f3]
H = straight_line_homotopy(F,G,t)
point = [CCi(1),CCi(-1), CCi(-1)]
A = jacobian_inverse(G,point)


    ## --- Initialization --------------------------------------------------
    coeff_ring  = base_ring(H[1])
    CCi         = coefficient_ring(coeff_ring)
    t           = 0
    h           = 1 / 10
    iter        = 0
    tprev       = 0
    n           = length(x)

    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
    A = jacobian_inverse(G, x)

    refinement_threshold = 1/8
    ## --- First Step (Linear predictor) ----------------------------------
    x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold)

    xprev  = x
    vprev  = v
    hprev  = h
    t     += h

    Ft, x, r, A, v, h, radii = refine_step(H, t, x, r, A, h; threshold = refinement_threshold)
    X  = hermite_predictor(H, x, xprev, v, vprev, hprev)
    tm = krawczyk_operator_taylor_model(H, X, t, A, r)


function truncate_degree(expr::T, var::String, max_degree::Int) where T <: AbstractAlgebra.Generic.MPoly
    terms = collect(AbstractAlgebra.terms(expr))
    vars = AbstractAlgebra.vars(expr)
    var_index = findfirst(x -> string(x) == var, vars)
    
    if isnothing(var_index)
        return expr
    end
    
    filtered_terms = filter(term -> 
        AbstractAlgebra.degree(term, var_index) <= max_degree, 
        terms)
    
    return sum(filtered_terms)
end    

    x,it = track(H, point, .1; iterations_count = true)