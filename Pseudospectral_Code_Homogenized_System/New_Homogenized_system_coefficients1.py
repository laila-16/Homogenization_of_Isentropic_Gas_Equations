from scipy import integrate
def double_bracket(f):
    # Compute the double bracket of a function
    mean = integrate.quad(f,0,1)[0]
    brace = lambda y: f(y)-mean
    brack_nzm = lambda y: integrate.quad(brace,0,y)[0]
    mean_bracket = integrate.quad(brack_nzm,0,1)[0]
    def brack(y):
        return integrate.quad(brace,0,y)[0] - mean_bracket
    return brack

def C_values(a, dady,p_0, delta):
    ainvsquared = lambda y: 1/a(y)**2
    ay_ainv = lambda y: ainvsquared(y)*dady(y)
    ainv = lambda y: 1/a(y)
    acubed=  lambda y: 1/a(y)**3
    ainvsquared = lambda y: 1/a(y)**2
    db_ay_ainv = double_bracket(ay_ainv)
    db_a = double_bracket(a)
    A= lambda y: a(y)*p_0
    Ay=lambda y: dady(y)*p_0
    A_inv1 = lambda y: 1/A(y)
    A_inv2 = lambda y: 1/A(y)**2
    A_inv3 = lambda y: 1/A(y)**3
    A_inv4 = lambda y: 1/A(y)**4
   #<a^-1[a]>
    db = double_bracket(a)
    integrand = lambda y: ainv(y)*db(y)
    C1 =integrate.quad(integrand,0,1)[0]

    #<a^-1 [[a]]>
    db = double_bracket(a)
    dbdb = double_bracket(db)
    integrand = lambda y: ainv(y)*dbdb(y)
    C2 =integrate.quad(integrand,0,1)[0]

    
    ay_aAinv2 =lambda y: dady(y)*ainv(y)*A_inv2(y)
    db =double_bracket(ay_aAinv2)
    integrand = lambda y: db(y)*a(y)
    C3=integrate.quad(integrand,0,1)[0]


    Ay_aAinv2 =lambda y: Ay(y)*ainv(y)*A_inv2(y)
    db =double_bracket(Ay_aAinv2)
    dbdb= double_bracket(db)
    integrand = lambda y: dbdb(y)*a(y)
    C4=integrate.quad(integrand,0,1)[0]


    Ay_aAinv2 =lambda y: Ay(y)*ainv(y)*A_inv2(y)
    db =double_bracket(Ay_aAinv2)
    integrand = lambda y: db(y)*A_inv1(y)
    C5=integrate.quad(integrand,0,1)[0]

    
    Ay_aAinv2 =lambda y: Ay(y)*ainv(y)*A_inv2(y)
    db =double_bracket(Ay_aAinv2)
    db1 = lambda y:db(y)*a(y)
    dbdb= double_bracket(db1)
    integrand = lambda y: dbdb(y)*Ay_aAinv2(y)
    C6=integrate.quad(integrand,0,1)[0]


    db = double_bracket(ainv)
    db1 = lambda y: db(y)*a(y)
    dbdb = double_bracket(db1)
    integrand = lambda y: dbdb(y)*Ay(y)*ainv(y)*A_inv2(y)
    C7=integrate.quad(integrand,0,1)[0]


    db = double_bracket(a)
    dba_Ay_aAinv2= lambda y: db(y)*Ay(y)*ainv(y)*A_inv2(y)
    dbdb= double_bracket(dba_Ay_aAinv2)
    integrand = lambda y: dbdb(y)*a(y)
    C8=integrate.quad(integrand,0,1)[0]


    db = double_bracket(a)
    dba_ainv= lambda y : db(y)*ainv(y)
    dbdb= double_bracket(dba_ainv)
    integrand = lambda y: dbdb(y)*a(y)
    C9=integrate.quad(integrand,0,1)[0]


    integrand = lambda y: a(y)*A_inv2(y)
    C10=integrate.quad(integrand,0,1)[0]


    db = double_bracket(ainv)
    db1 = lambda y: db(y)*a(y)
    dbdb = double_bracket(db1)
    integrand = lambda y: dbdb(y)*ainv(y)
    C11=integrate.quad(integrand,0,1)[0]


    integrand = lambda y: ainv(y)*A_inv2(y)
    C12=integrate.quad(integrand,0,1)[0]


    Ay_aAinv2 =lambda y: Ay(y)*ainv(y)*A_inv2(y)
    db =double_bracket(Ay_aAinv2)
    integrand = lambda y: db(y)*a(y)
    C13=integrate.quad(integrand,0,1)[0]


    integrand = lambda y: A_inv1(y)
    C14=integrate.quad(integrand,0,1)[0]


    db = double_bracket(A_inv1)
    integrand = lambda y: db(y)*a(y)
    C15=integrate.quad(integrand,0,1)[0]


    #<a>^2
    a2 = lambda y: a(y)**2
    avga=(1/delta)*integrate.quad(a,0,1)[0]
    avga2= (1/delta)*integrate.quad(a2,0,1)[0]
    ainvavg = (1/delta)*integrate.quad(ainv,0,1)[0]
    return C1, C2, C3, C4, C5, C6, C7, C8, C9,C10,C11,C12,C13,C14,C15,avga,ainvavg
    
    
def Homogenized_system_coef(C1, C2, C3, C4, C5, C6, C7, C8, C9,C10,C11,C12,C13,C14,C15,avga,ainvavg, p_0, P1, P11, delta):
    r1 = - 1/avga
    r2 = delta* (-C1/(ainvavg* avga**2))
    r3 = delta* (2*C13/(ainvavg*avga))
    r4 = delta**2 *(-C3/(P1*avga**2) + (4*C13**2)/(ainvavg*P1*avga**2) +(4*C13*C14)/(ainvavg*P1*avga**2) - (C13*P11)/(P1**2 * avga**2))
    r5 = delta**2 *((C9)/(ainvavg*avga**3) - C2/(ainvavg*avga**2))
    r5b = delta**2 *((C9)/(ainvavg*avga**3) - C2/(ainvavg*avga**2))*(-1/avga)
    r6 = delta**2 * (-2*C3/(ainvavg*avga))
    r7 = delta**2 *(((2*C8)/(ainvavg*avga**2)) -(2*C4)/(ainvavg*avga))
    r8 = delta**2 * (((-2*C1*C13)/(ainvavg**2 * avga**2))+2*C8/(ainvavg*avga**2) -2*C4/(ainvavg*avga) )
    beta1 = -P1 /ainvavg
    beta2 = delta* 2*(-C13-C14)/(ainvavg*avga) 
    beta3 = -delta* P11/ainvavg
    beta4 = delta* C1*P1/(avga*ainvavg**2)
    beta5 = delta* (2*C3+2*C10)/(ainvavg*avga)
    beta6 = delta**2 *((2*C1*C13 )/(ainvavg**2 * avga**2)+(-C8)/(ainvavg*avga**2))
    beta7 = delta**2 * (-(2*C1*C13)/(ainvavg**2 * avga**2)-(2*C4/ainvavg*avga)+(4*C7/ainvavg**2 *avga))
    beta8 = delta**2 *((C12-3*C5+4*C6)/(ainvavg**2)+(4*C13**2 +4*C13*C14)/(ainvavg**2 *avga))
    beta9 = delta**2 *((2*C7*P1/ainvavg**3)-(2*C1*C13*P1)/(avga*ainvavg**3)+(C1*P11)/(avga*ainvavg**2))
    beta10 = delta**2 *(C1*P11/(avga*ainvavg**2))
    beta11= delta**2 *(C11*P1/(avga*ainvavg**3) - C2*P1/(avga*ainvavg**2))
    beta11b = delta**2 *(C11*P1/(avga*ainvavg**3) - C2*P1/(avga*ainvavg**2))*(-P1/ainvavg)
    return r1, r2, r3, r4, r5,r5b, r6,r7,r8, beta1, beta2, beta3, beta4, beta5, beta6, beta7,beta8, beta9, beta10, beta11,beta11b







