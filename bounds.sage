import itertools

def mathem(ineq):
    r = "{}".format(ineq)
    z = ""
    for s in r:
        s = s.replace('T','').replace('\r','').replace('<','').replace('=','')
        z+=s
    return(float(sage_eval(z)))
print mathem('T <= 1/2*sqrt(10) - 1/4*sqrt(2)')

angles = [pi/2, pi/3, pi/4, pi/6]

B = itertools.product(*[angles for i in range(5)])

A.<T> = PolynomialRing(QQ)

sols = []
CH = {}


for b in B:
    if ((b[0]+b[1]+b[2] - pi > 0) and (b[0]+b[3]+b[4] - pi > 0)):
        #G = matrix([[1, -cos(b[0]), -cos(b[1]), -cos(b[3])], [-cos(b[0]), 1, -cos(b[2]), -cos(b[4])], [-cos(b[1]), -cos(b[2]), 1, -T], [-cos(b[3]), -cos(b[4]), -T, 1]])
        G1 = matrix([[1, -cos(b[0]), -cos(b[1])], [-cos(b[0]), 1, -cos(b[2])], [-cos(b[3]), -cos(b[4]), -T]])
        G2 = matrix([[1, -cos(b[0]), -cos(b[3])], [-cos(b[0]), 1, -cos(b[4])], [-cos(b[3]), -cos(b[4]), 1]])
        G3 = matrix([[1, -cos(b[0]), -cos(b[1])], [-cos(b[0]), 1, -cos(b[2])], [-cos(b[1]), -cos(b[2]), 1]])
        
        a1 = (cos(b[1]) + cos(b[0])*cos(b[2]))/(sin(b[0])*sin(b[2]))
        a3 = (cos(b[2]) + cos(b[0])*cos(b[1]))/(sin(b[0])*sin(b[1]))
        
        a2 = (cos(b[3]) + cos(b[0])*cos(b[4]))/(sin(b[0])*sin(b[4]))
        a4 = (cos(b[4]) + cos(b[0])*cos(b[3]))/(sin(b[0])*sin(b[3]))
        
        a1s = sqrt((1-a1)/2)
        a1c = sqrt((1+a1)/2)
        
        a2s = sqrt((1-a2)/2)
        if (a2s == 0):
            print b
            print b[0]+b[3]+b[4]
        
        a2c = sqrt((1+a2)/2)
        
        a3s = sqrt((1-a3)/2)
        a3c = sqrt((1+a3)/2)
        
        a4s = sqrt((1-a4)/2)
        a4c = sqrt((1+a4)/2)
        
        A0 = tanh(ln(cot(b[0]/4)))

        lng1 = arcsinh(A0/(a1s/a1c)) + arcsinh(A0/(a2s/a2c))
        lng2 = arcsinh(A0/(a3s/a3c)) + arcsinh(A0/(a4s/a4c))
        
        l = max(lng1,lng2)
        
        T = var('T')
        eq = (- G1.det(A)/sqrt(abs(det(G2)*det(G3))) <= float(cosh(l)))
        #sols.append(solve(eq, T, solution_dict=True)[0])
        
        st = solve(eq, T, solution_dict=True)[0][0]
        #print st
        CH[(float(b[0]), float(b[1]), float(b[2]), float(b[3]), float(b[4]))] = mathem(st)
        sols.append(mathem(st))
#float(b[0]), float(b[1]), float(b[2]), float(b[3]), float(b[4])
print CH[(float(pi/2), float(pi/2), float(pi/2), float(pi/2), float(pi/2))]
