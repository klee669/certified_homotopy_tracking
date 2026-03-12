using Pkg
Pkg.activate(joinpath(@__DIR__, "../..")) 

using Nemo
using AbstractAlgebra
using CertifiedMonodromyComputation
using GAP

println("=== Running Mathieu Group Example (Manual Loops) ===")

# ------------------------------------------------------------------------------
# 1. Setup System
# ------------------------------------------------------------------------------
@monodromy_setup begin
    vars = (x)
    params = (h)
end
const CCi = _CCi


g0 = CCi(.549472304541888876273331247256+.675650335767895613928004918526*im)
g = g0
tau =  (3^17) * (47323*g0^3 - 1084897*g0^2 + 7751*g0 - 711002) / (23^3)
tau = (2^38) * tau
minPoly = g^4 + g^3 + 9*g^2 - 10*g + 8
P2 = (8*g^3 + 16 * g^2 - 20*g + 20) * x^2 - (7*g^3 + 17*g^2 - 7*g + 76) * x  + (-13*g^3 + 25*g^2 - 107*g + 596)
P3  = 8*(31*g^3 + 405*g^2 - 459*g + 333) * x^3 + (941*g^3 + 1303 * g^2 -1853*g + 1772)*x +(85*g^3 - 385*g^2 + 395*g - 220)
P4 = 32 * (4*g^3 - 69 * g^2 + 74*g - 49) * x^4 + 32*(21*g^3 + 53*g^2 - 68*g + 58)*x^3 - 8*(97*g^3+95*g^2-145*g+148)*x^2 + 8 * (41*g^3 - 89*g^2-g + 140)*x +(-123*g^3+391*g^2-93*g+3228)
P42= P4^2
P44= P42^2
P22 = P2^2
fBelyi = P22*P44*P3;
fBelyi/tau-1


f = [fBelyi/tau-h];
bp = [CCi(1/2)]
# ------------------------------------------------------------------------------
# 2. Local Helper Functions
# ------------------------------------------------------------------------------
function search_point_local(res, p_list)
    n = length(p_list);
    k = 1;
    dummy = CertifiedMonodromyComputation.max_norm(matrix(res-p_list[1]));
    for i = 2:n 
        m = CertifiedMonodromyComputation.max_norm(matrix(res-p_list[i]));
        if m < dummy
            dummy = m;
            k = i;
        end
    end
    k
end
function track_loop(bp, a, b, x0, r, p_list, i, F)
    println("Root Number $i: Tracking the first edge")
    F1 = specified_system(bp, a, F);
    x1 = track(F1,x0; show_display=false);
    if x1 === nothing return nothing, nothing end 

    println("Root Number $i: Tracking the second edge")
    F2 = specified_system(a, b, F);
    x2 = track(F2,x1; show_display=false);
    if x2 === nothing return nothing, nothing end 

    println("Root Number $i: Tracking the third edge")
    F3 = specified_system(b, bp, F);
    x3 = track(F3,x2; show_display=false);
    if x3 === nothing return nothing, nothing end 

    ind = search_point_local(x3, p_list);
    println("Result: Mapped to $ind")
    return x3, ind
end

function generate_perm(F, bp, a, b, r, p_list)
    n = length(p_list);
    perm = [];
    for i = 1:n
        res, ind = track_loop(bp,a,b,p_list[i],r, p_list, i, F);
        
        if res === nothing 
            println("Stopped by user. Returning partial permutation.")
            break 
        end
        
        push!(perm, ind);
    end
    perm
end

# ------------------------------------------------------------------------------
# 3. Manual Loop Definitions
# ------------------------------------------------------------------------------
r = 0.1
p_list=[[CCi(-8.46348528726219075641182680047306619070122436150937232942857413267516074415,-8.236051504540776298957760979865031236775465557582439518026284354264711900635)]]
# red loop
red1 = [CCi(1/2)]
red2 = [CCi(-1,1)]
red3 = [CCi(-1,-1)]
println("\n[Loop 1] Red Loop")
p1 = generate_perm(f, red1, red2, red3, r, p_list) #[8, 2, 3, 4, 24, 6, 7, 1, 27, 23, 11, 16, 13, 14, 25, 12, 17, 18, 19, 20, 21, 22, 10, 5, 15, 26, 9]


# green loop
green1 = [-500,0]
green2 = [-500,CCi(-800,-500)]
green3 = [-500,CCi(-800,500)]
println("\n[Loop 2] Green Loop")
p2 = generate_perm(f, green1, green2, green3, r, p_list) #[25, 4, 17, 2, 27, 11, 26, 15, 24, 12, 6, 10, 14, 13, 8, 23, 3, 18, 19, 22, 21, 20, 16, 9, 1, 7, 5]


# green+purple loop
gp1 = [-500,0]
gp2 = [-500,CCi(-1320,-500)]
gp3 = [-500,CCi(-1320,500)]
println("\n[Loop 3] Green+Purple Loop")
p3 = generate_perm(f, gp1, gp2, gp3, r, p_list) #[8, 2, 3, 4, 24, 6, 7, 1, 27, 23, 11, 16, 13, 14, 25, 12, 17, 18, 19, 20, 21, 22, 10, 5, 15, 26, 9]

# green+purple+blue loop
gpb1 = [-500,0]
gpb2 = [-500,CCi(-3000,-500)]
gpb3 = [-500,CCi(-3000,500)]
println("\n[Loop 4] Green+Purple+Blue Loop")
p4 = generate_perm(f, gpb1, gpb2, gpb3, r, p_list) #[8, 2, 3, 4, 24, 6, 7, 1, 27, 23, 11, 16, 13, 14, 25, 12, 17, 18, 19, 20, 21, 22, 10, 5, 15, 26, 9]


# ------------------------------------------------------------------------------
# 4. GAP Analysis
# ------------------------------------------------------------------------------
println("\n=== GAP Analysis ===")

p1_gap = GAP.Globals.PermList(GAP.Obj(p1))
p2_gap = GAP.Globals.PermList(GAP.Obj(p2))

prod = p1_gap * p2_gap
println("p1 * p2 calculated in GAP.")

G = GAP.Globals.Group(p1_gap, p2_gap)
println("Group G defined.")

println("Structure Description:")
println(GAP.Globals.StructureDescription(G)) # C2 x C2 (= K4)

println("Galois Width:")
gw = galois_width(G) # 2
