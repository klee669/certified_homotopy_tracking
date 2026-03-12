include("../src/certified_monodromy_computation.jl")

# polynomial rings
@setupfield begin
    AcbField()
    (x,y)
    (η,)
    (t,)
    (a,)
end
CCi = _CCi



# system
f1= 121/21+13/11*x^4-15/13*x^3*y^2-(13/11)*x^2*y^3-11/14*y^4+27/21*x^2*y^2
f2= 14/5+33/15*x-27/13*y+a*x^2*y^2

f = [_PR(f1) _PR(f2)]
p_list = toCCi([[-.374778+.679561*im, .02235+2.10128*im], [.553921-1.14395*im, .538205-1.04926*im], [-1.12533+.876307*im, -.259731+.593236*im], [-1.04353-.200625*im, -.294563-1.25303*im], [1.06328-1.00602*im,.43207+1.5035*im], [1.41768-.006179*im, -.622983-1.63078*im], [.216215-.0486911*im, 1.63823+.003928*im], [.501993-1.07923*im, -.074731+2.05441*im], [-.820392-.472532*im, -1.62938-1.1245*im],[-.394116+.649671*im, 1.55253+.075212*im], [1.23719+1.34148*im, .901466+.402221*im], [-.941442+.861031*im, 1.16403-.440796*im], [1.08427+.741466*im, -1.70977-1.17152*im], [-1.08474-1.69364*im,-.834192+.567681*im], [-1.06465-1.27026*im, .479207-.521953*im], [.774428+1.77162*im, -1.30274-.10964*im]],CCi)

function search_point(res, p_list)
    n = length(p_list);
    k = 1;
    dummy = max_norm(matrix(res-p_list[1]));
    for i = 2:n 
        m = max_norm(matrix(res-p_list[i]));
        if m < dummy
            dummy = m;
            k = i;
        end
    end
    k
end

function track_loop(bp, a, b, c, d, e, x0, r, p_list, i, F)
    println("Root Number $i: Tracking the first edge")
    F1 = specified_system(bp, a, F);
    x1 = track(F1,x0,r);

    println("Root Number $i: Tracking the second edge")
    F2 = specified_system(a, b, F);
    x2 = track(F2,x1,r);

    println("Root Number $i: Tracking the third edge")
    F3 = specified_system(b, c, F);
    x3 = track(F3,x2,r);

    println("Root Number $i: Tracking the third edge")
    F4 = specified_system(c, d, F);
    x4 = track(F4,x3,r);

        println("Root Number $i: Tracking the third edge")
    F5 = specified_system(d, e, F);
    x5 = track(F5,x4,r);

        println("Root Number $i: Tracking the third edge")
    F6 = specified_system(e, bp, F);
    x6 = track(F6,x5,r);

    println( search_point(x6, p_list));
    x6, search_point(x6, p_list)
end

function generate_perm(F, bp, a, b, c, d, e, r, p_list)
    n = length(p_list);
    perm = [];
    for i = 1:n
        res, ind = track_loop(bp, a, b, c, d, e,p_list[i],r, p_list, i, F);
        perm = push!(perm, ind);
    end
    perm
end



# red loop
point = [CCi(.5,1.2)]
point1 = [CCi(.01,.01)]
point2 = [CCi(.01,-.01)]
point3 = [CCi(-.01,-.01)]
point4 = [CCi(-.01,.01)]
point5 = [CCi(.01,.01)]

r = 0.1
p1 = generate_perm(f, point, point1, point2, point3, point4, point5, r, p_list) #[13, 2, 3, 14, 5, 8, 1, 10, 7, 11, 6, 12, 4, 16, 15, 9]
p1 = string(p1)

using GAP
p1_str = "p1:= PermList([9, 2, 3, 4, 13, 12, 7, 6, 5, 10, 14, 8, 1, 16, 11, 15])"

GAP.evalstr(p1_str) #p1 is a product of 12 transpositions






# system2
f1= 121/21+13/11*x^4-151/113*x^3*y^2-(13/11)*x^2*y^3-15/14*y^4+27/21*x^2*y^2
f2= 14/25+33/15*x-127/131*y+a*x^2*y^2

f = [_PR(f1) _PR(f2)]
p_list = toCCi([[-.489935-.33349*im, -1.54386-.217135*im], [-.156376+.417507*im, .025611+1.6079*im], [-.496276-.266691*im, -.141091-1.54464*im], [1.03419-.942003*im, .38555+1.2045*im], [1.14014+1.37107*im,.79025+.458786*im], [.924465+.954986*im, -1.26061-.67921*im], [-.927714+.686369*im, 1.3477-.332334*im], [.332293-.953761*im, -.057588+1.73955*im], [.572033-.826477*im, .198786-1.14487*im],[.30581-.0806357*im, 1.51187+.007455*im], [1.09905-.233895*im, -.364091-1.40025*im], [-.398706+.930409*im, 1.34213+.20594*im], [.501317+1.37852*im, -1.22197-.231167*im], [-1.37098-1.4837*im,-.597223+.680622*im], [-1.27782-1.42898*im, .422884-.756778*im], [-.791493+.810779*im, -.83835+.401633*im]],CCi)


# red loop
point = [CCi(.5,1.2)]
point1 = [CCi(.01,.01)]
point2 = [CCi(.01,-.01)]
point3 = [CCi(-.01,-.01)]
point4 = [CCi(-.01,.01)]
point5 = [CCi(.01,.01)]

r = 0.1
p1 = generate_perm(f, point, point1, point2, point3, point4, point5, r, p_list) #[13, 2, 3, 14, 5, 8, 1, 10, 7, 11, 6, 12, 4, 16, 15, 9]
p1 = string(p1)

using GAP
p1_str = "p1:= PermList([1, 2, 4, 11, 14, 7, 8, 6, 9, 10, 12, 3, 15, 13, 5, 16])"

GAP.evalstr(p1_str) #p1 is a product of 12 transpositions


# system3
f1= -12/21+13/11*x^4-151/113*x^3*y^2+13/11*x^2*y^3+15/14*y^4-27/21*x^2*y^2
f2= 14/25+33/15*x-127/131*y+a*x^2*y^2

f = [_PR(f1) _PR(f2)]
p_list = toCCi([[.119821-.0055793*im, .856181-.000170*im], [-.909305+.007889*im, -.773233+.312059*im], [-1.32591-.306618*im, .758219-.843046*im], [-.482307+.794729*im, -.941802+.524969*im], [1.19071+.404223*im,.824311+.695191*im], [-.378375-.321542*im, -.07274-.910939*im], [1.41643+1.51054*im, -.898406-.507628*im], [.552007+1.13135*im, .935431+.453169*im], [.51531-.864817*im, .21679-1.11382*im],[.429433-.837128*im, .182548+1.81777*im], [-.210194+.342144*im, -.020766+.906012*im], [-.562096+.850074*im, 1.42555+.052282*im], [.845681+.452309*im, -1.38504-1.3276*im], [-.406215-1.87971*im,-.414854+.936459*im], [-.212757-1.04221*im, .311543-.927854*im], [-.582231-.235649*im, -1.00373-.066846*im]],CCi)


# red loop
point = [CCi(.5,1.2)]
point1 = [CCi(.001,.001)]
point2 = [CCi(.001,-.001)]
point3 = [CCi(-.001,-.001)]
point4 = [CCi(-.001,.001)]
point5 = [CCi(.001,.001)]

r = 0.1
p1 = generate_perm(f, point, point1, point2, point3, point4, point5, r, p_list) #[1, 2, 8, 9, 4, 6, 3, 14, 5, 13, 11, 15, 12, 7, 10, 16]
p1 = string(p1)

using GAP
p1_str = "p1:= PermList([1, 2, 8, 9, 4, 6, 3, 14, 5, 13, 11, 15, 12, 7, 10, 16])"

GAP.evalstr(p1_str) #p1 is a product of 12 transpositions

