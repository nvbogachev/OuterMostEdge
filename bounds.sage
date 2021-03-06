import itertools

def mathem(ineq):
    r = "{}".format(ineq)
    z = ""
    for s in r:
        s = s.replace('T','').replace('\r','').replace('<','').replace('=','')
        z+=s
    return(float(sage_eval(z)))
print mathem('T <= 1/2*sqrt(10) - 1/4*sqrt(2)')

import itertools
from collections import OrderedDict

#
# Possible angles for proof of Main Theorem 1 after Propositions 3.1 and 3.2
#

angles = [pi/2, pi/3, pi/4, pi/5]

B = itertools.product(*[angles for i in range(5)])

A.<T> = PolynomialRing(QQ)

sols = []
CH = {}
bounds = {}


# b[0] = alpha_12
# b[1] = alpha_13
# b[2] = alpha_23
# b[3] = alpha_14
# b[4] = alpha_24

for b in B:
    #
    # Check Lemma 2.1 (i) condition:
    #
    if ((b[0]+b[1]+b[2] - pi > 0) and (b[0]+b[3]+b[4] - pi > 0) and b[0]<pi/2): ### b[0] = angle_12 < pi/2

        #
        # Gram matrix G = G(u1,u2,u3,u4)
        #
        # G = matrix([[1, -cos(b[0]), -cos(b[1]), -cos(b[3])], [-cos(b[0]), 1, -cos(b[2]), -cos(b[4])], [-cos(b[1]), -cos(b[2]), 1, -T], [-cos(b[3]), -cos(b[4]), -T, 1]])
        #
        # Algebraic complements G_ij:
        #

        G34 = matrix([[1, -cos(b[0]), -cos(b[1])], [-cos(b[0]), 1, -cos(b[2])], [-cos(b[3]), -cos(b[4]), -T]])
        G33 = matrix([[1, -cos(b[0]), -cos(b[3])], [-cos(b[0]), 1, -cos(b[4])], [-cos(b[3]), -cos(b[4]), 1]])
        G44 = matrix([[1, -cos(b[0]), -cos(b[1])], [-cos(b[0]), 1, -cos(b[2])], [-cos(b[1]), -cos(b[2]), 1]])

        a1 = (cos(b[1]) + cos(b[0])*cos(b[2]))/(sin(b[0])*sin(b[2])) # cos alpha_1
        a3 = (cos(b[2]) + cos(b[0])*cos(b[1]))/(sin(b[0])*sin(b[1])) # cos alpha_3

        a2 = (cos(b[3]) + cos(b[0])*cos(b[4]))/(sin(b[0])*sin(b[4])) # cos alpha_2
        a4 = (cos(b[4]) + cos(b[0])*cos(b[3]))/(sin(b[0])*sin(b[3])) # cos alpha_4

        a1s = sqrt((1-a1)/2) # sin (alpha_1/2)
        a1c = sqrt((1+a1)/2) # cos (alpha_1/2)

        a2s = sqrt((1-a2)/2) # sin (alpha_2/2)
        a2c = sqrt((1+a2)/2) # cos (alpha_2/2)

        a3s = sqrt((1-a3)/2) # sin (alpha_3/2)
        a3c = sqrt((1+a3)/2) # cos (alpha_3/2)
        
        a4s = sqrt((1-a4)/2) # sin (alpha_4/2)
        a4c = sqrt((1+a4)/2) # cos (alpha_4/2)

        A0 = tanh(ln(cot(b[0]/4))) # this is actually cos (b[0]/2) = cos (alpha_12/2)
        
        lng1 = arcsinh(A0/(a1s/a1c)) + arcsinh(A0/(a2s/a2c)) # F_{1,2} (alpha)
        lng2 = arcsinh(A0/(a3s/a3c)) + arcsinh(A0/(a4s/a4c)) # F_{3,4} (alpha)

        l = max(lng1,lng2) # Corollary 3: F(alpha) = max {F_{1,2} (alpha), F_{3,4} (alpha)}
        
        T = var('T') # width of a small ridge
        #
        # Inequality: G34/sqrt(G33*G44) = cosh(a) <= cosh F(alpha)
        #
        eq = (- G34.det(A)/sqrt(abs(det(G33)*det(G44))) <= float(cosh(l)))
        #sols.append(solve(eq, T, solution_dict=True)[0])
        
        st = solve(eq, T, solution_dict=True)[0][0]

        # CH[(float(b[0]), float(b[1]), float(b[2]), float(b[3]), float(b[4]))] = mathem(st)
        bounds[(b[0],b[1],b[2],b[3],b[4])] = mathem(st)
        sols.append(mathem(st))
# (b[0],b[1],b[2],b[3],b[4])
#print 'CH=', CH[(float(pi/2), float(pi/2), float(pi/2), float(pi/2), float(pi/2))]

print bounds

#print(latex(bounds))

#for t in sorted(sols):
#    print t

bbbb = len(list(OrderedDict.fromkeys(sols)))
print 'Number of different ridges =', bbbb
