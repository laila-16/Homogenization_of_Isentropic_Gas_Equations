import numpy as np
from scipy import integrate

def double_bracket(f):
    # Compute the mean of the function over [0, 1]
    mean_1 = integrate.quad(f, 0, 0.5, limit=100)[0]
    mean_2 = integrate.quad(f, 0.5, 1, limit=100)[0]
    mean = mean_1 + mean_2
    
    # Define the brace function
    brace = lambda y: f(y) - mean
    
    # Define the non-zero mean bracket function
    def brack_nzm(y):
        if y <= 0.5:
            return integrate.quad(brace, 0, y, limit=100)[0]
        else:
            integral_1 = integrate.quad(brace, 0, 0.5, limit=100)[0]
            integral_2 = integrate.quad(brace, 0.5, y, limit=100)[0]
            return integral_1 + integral_2
    
    # Compute the mean of the bracket function over [0, 1]
    mean_bracket_1 = integrate.quad(brack_nzm, 0, 0.5, limit=100)[0]
    mean_bracket_2 = integrate.quad(brack_nzm, 0.5, 1, limit=100)[0]
    mean_bracket = mean_bracket_1 + mean_bracket_2
    
    # Define the final bracket function
    def brack(y):
        if y <= 0.5:
            return integrate.quad(brace, 0, y, limit=100)[0] - mean_bracket
        else:
            integral_1 = integrate.quad(brace, 0, 0.5, limit=100)[0]
            integral_2 = integrate.quad(brace, 0.5, y, limit=100)[0]
            return integral_1 + integral_2 - mean_bracket
    
    return brack

def C_values(a, dady, delta):
    ainvsquared = lambda y: 1/a(y)**2
    ay_ainv = lambda y: ainvsquared(y) * dady(y)
    ainv = lambda y: 1/a(y)
    acubed = lambda y: 1/a(y)**3
    
    # Compute double brackets
    db_ay_ainv = double_bracket(ay_ainv)
    db_a = double_bracket(a)
    
    # Compute various averages
    inv_avg_1 = integrate.quad(ainv, 0, 0.5, limit=100)[0]
    inv_avg_2 = integrate.quad(ainv, 0.5, 1, limit=100)[0]
    inv_avg = inv_avg_1 + inv_avg_2
    
    a_avg_1 = integrate.quad(a, 0, 0.5, limit=100)[0]
    a_avg_2 = integrate.quad(a, 0.5, 1, limit=100)[0]
    a_avg = a_avg_1 + a_avg_2
    
    a_invsq_1 = integrate.quad(ainvsquared, 0, 0.5, limit=100)[0]
    a_invsq_2 = integrate.quad(ainvsquared, 0.5, 1, limit=100)[0]
    a_invsq = a_invsq_1 + a_invsq_2
    
    acubed_avg_1 = integrate.quad(acubed, 0, 0.5, limit=100)[0]
    acubed_avg_2 = integrate.quad(acubed, 0.5, 1, limit=100)[0]
    acubed_avg = acubed_avg_1 + acubed_avg_2
    
    p_0 = 0.3
    A = lambda y: a(y) * p_0
    Ay = lambda y: dady(y) * p_0
    A_inv1 = lambda y: 1/A(y)
    A_inv2 = lambda y: 1/A(y)**2
    A_inv3 = lambda y: 1/A(y)**3
    A_inv4 = lambda y: 1/A(y)**4
    
    # Compute coefficients
    db = double_bracket(a)
    integrand = lambda y: ainv(y) * db(y)
    C1_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C1_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C1 = C1_1 + C1_2

    db = double_bracket(a)
    dbdb = double_bracket(db)
    integrand = lambda y: ainv(y) * dbdb(y)
    C2_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C2_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C2 = C2_1 + C2_2

    C3 = -0.5 * (1 / p_0**2) * (inv_avg - a_avg * a_invsq)

    db = double_bracket(ainvsquared)
    integrand = lambda y: db(y) * a(y)
    C4_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C4_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C4 = -0.5 * (1 / p_0) * (C4_1 + C4_2)

    C5 = -0.5 * (1 / p_0**2) * (acubed_avg - inv_avg * a_invsq)

    C6 = -0.25 * (1 / p_0**2) * (acubed_avg + a_avg * a_invsq**2 - 2 * inv_avg * a_invsq)

    C7 = 0.5 * (1 / p_0) * a_invsq * C1

    C8 = 0

    db = double_bracket(a)
    dba_ainv = lambda y: db(y) * ainv(y)
    dbdb = double_bracket(dba_ainv)
    integrand = lambda y: dbdb(y) * a(y)
    C9_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C9_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C9 = C9_1 + C9_2

    integrand = lambda y: ainv(y)
    C10_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C10_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C10 = C10_1 + C10_2

    db = double_bracket(ainv)
    db1 = lambda y: db(y) * a(y)
    dbdb = double_bracket(db1)
    integrand = lambda y: dbdb(y) * ainv(y)
    C11_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C11_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C11 = C11_1 + C11_2

    integrand = lambda y: ainv(y) * A_inv2(y)
    C12_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C12_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C12 = C12_1 + C12_2

    Ay_aAinv2 = lambda y: Ay(y) * ainv(y) * A_inv2(y)
    db = double_bracket(Ay_aAinv2)
    integrand = lambda y: db(y) * a(y)
    C13_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C13_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C13 = C13_1 + C13_2

    integrand = lambda y: A_inv1(y)
    C14_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C14_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C14 = C14_1 + C14_2

    db = double_bracket(A_inv1)
    integrand = lambda y: db(y) * a(y)
    C15_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C15_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C15 = C15_1 + C15_2

    integrand = lambda y: a(y) * A_inv2(y)
    C16_1 = integrate.quad(integrand, 0, 0.5, limit=100)[0]
    C16_2 = integrate.quad(integrand, 0.5, 1, limit=100)[0]
    C16 = C16_1 + C16_2

    # Compute <a> and <a^-1>
    a2 = lambda y: a(y)**2
    avga_1 = (1/delta) * integrate.quad(a, 0, 0.5, limit=100)[0]
    avga_2 = (1/delta) * integrate.quad(a, 0.5, 1, limit=100)[0]
    avga = avga_1 + avga_2
    
    avga2_1 = (1/delta) * integrate.quad(a2, 0, 0.5, limit=100)[0]
    avga2_2 = (1/delta) * integrate.quad(a2, 0.5, 1, limit=100)[0]
    avga2 = avga2_1 + avga2_2
    
    ainvavg_1 = (1/delta) * integrate.quad(ainv, 0, 0.5, limit=100)[0]
    ainvavg_2 = (1/delta) * integrate.quad(ainv, 0.5, 1, limit=100)[0]
    ainvavg = ainvavg_1 + ainvavg_2
    
    return C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16, avga, ainvavg


def Homogenized_system_coef(C1, C2, C3, C4, C5, C6, C7, C8, C9,C10,C11,C12,C13,C14,C15,C16,avga,ainvavg, p_0, P1, P11, delta):
    r1 = - 1/avga
    r2 = delta* (-C1/(C10* avga**2))
    r3 = delta* (2*C13/(C10*avga))
    r4 = delta**2 *(-C3/(P1*avga**2) + (4*C13**2)/(C10*P1*avga**2) +(4*C13*C14)/(C10*P1*avga**2) - (C13*P11)/(P1**2 * avga**2))
    r5 = delta**2 *((C9)/(C10*avga**3) - C2/(C10*avga**2))
    r5b = delta**2 *((C9)/(C10*avga**3) - C2/(C10*avga**2))*(-1/avga)
    r6 = delta**2 * (-2*C3/(C10*avga))
    r7 = delta**2 *(((-2*C15+2*C8)/(C10*avga**2)) -(2*C1*C14)/(C10**2 * avga**2)-(2*C4)/(C10*avga))
    r8 = delta**2 * (((-2*C1*C13)/(C10**2 * avga**2))- 2*C1*C14/(C10**2 * avga**2) -2*C15/(C10*avga**2) +2*C8/(C10*avga**2) -2*C4/(C10*avga) )
    beta1 = -P1 /C10
    beta2 = delta* 2*(-C13-C14)/(C10*avga) 
    beta3 = -delta* P11/C10
    beta4 = delta* C1*P1/(avga*C10**2)
    beta5 = delta* (2*C3+2*C16)/(C10*avga)
    beta6 = delta**2 *((2*C1*C13 + 2*C1*C14)/(C10**2 * avga**2)+(2*C15-C8)/(C10*avga**2))
    beta7 = delta**2 * (-(2*C1*C13)/(C10**2 * avga**2)-(2*C4/C10*avga)+(4*C7/C10**2 *avga))
    beta8 = delta**2 *((C12-3*C5+4*C6)/(C10**2)+(4*C13**2 +4*C13*C14)/(C10**2 *avga))
    beta9 = delta**2 *((2*C7*P1/C10**3)-(2*C1*C13*P1)/(avga*C10**3)+(C1*P11)/(avga*C10**2))
    beta10 = delta**2 *(C1*P11/(avga*C10**2))
    beta11= delta**2 *(C11*P1/(avga*C10**3) - C2*P1/(avga*C10**2))
    beta11b = delta**2 *(C11*P1/(avga*C10**3) - C2*P1/(avga*C10**2))*(-P1/ainvavg)
    return r1, r2, r3, r4, r5,r5b, r6,r7,r8, beta1, beta2, beta3, beta4, beta5, beta6, beta7,beta8, beta9, beta10, beta11,beta11b



