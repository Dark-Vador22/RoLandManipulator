import sympy as sp
from sympy import symbols, pprint

N15, N56, N78, N27, N39, N911, N410, N1012 = symbols("N15, N56, N78, N27, N39, N911, N410, N1012")

g1 = N15*N56
g2 = N27*N78
g3 = (N56*N15*N78*N27 - N56*N15*N39)*N911
g4 = N78*N27*N911
g5 = N39*N911
g6 = (N56*N15*N78*N27 - N56*N15*N410)*N1012
g7 = N78*N27*N1012
g8 = N410*N1012

e1 = sp.simplify(g1*g4 - g3)/(g1*g5)
e2 = sp.simplify(g1*g7 - g6)/(g1*g8)
e3 = g4/(g2*g5)
e4 = g7/(g2*g8)

pprint(e1)
pprint(e2)
pprint(e3)
pprint(e4)

