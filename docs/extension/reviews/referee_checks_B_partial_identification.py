"""Independent checks. Run with python3 -I referee_checks.py (NumPy and SciPy)."""
from fractions import Fraction as F
from itertools import product

import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import linprog
from scipy.stats import norm

def cov(w, a, b):
    w = w / w.sum()
    return (w * (a - w @ a)) @ (b - w @ b)

def within(x, m1, m0):
    p = np.ones(len(x)) / len(x)
    mu = p @ x
    return mu*cov(p*x, x, m1) + (1-mu)*cov(p*(1-x), x, m0)

# M1. Direct weights with absent cells, including a pure first row.
X = np.array([[1,0,0],[0,.6,.4],[.3,0,.7],[.2,.5,.3]])
p = np.array([1,2,3,4])/10
m = np.array([[.8,.3,.1],[.2,.9,.5],[.3,.1,.6],[.1,.5,.9]])
T, mu = (X*m).sum(axis=1), p @ X
B = (p @ (X*m))/mu
M = X.T @ (p[:, None]*X)
W2 = sum(mu[s]*np.array([cov(p*X[:,s], X[:,r], m[:,s])
                         for r in range(3)]) for s in range(3))
assert np.allclose(X.T @ (p*T)-M@B, W2) and abs(W2.sum())<1e-12
for r in range(3):
    x = X[:,r]
    d = B[r] - (p@T-mu[r]*B[r])/(1-mu[r])
    other = np.divide(T-x*m[:,r], 1-x, out=np.zeros(len(x)), where=x!=1)
    W1 = mu[r]*cov(p*x,x,m[:,r])+(1-mu[r])*cov(p*(1-x),x,other)
    between = p @ (x*(m[:,r]-T))
    assert np.isclose(W1, cov(p,x,T)-cov(p,x,x)*d)
    assert np.isclose(between, mu[r]*(1-mu[r])*d-cov(p,x,T))
print("M1: all boundary-design identities hold")

# M2/M3. Enumerate individual counts for the stated three-unit example.
rows, yes = [1,2,3], [0,1,0]
legal = [range(max(0,y-(4-r)),min(r,y)+1) for r,y in zip(rows,yes)]
gaps = {F(sum(z),6)-F(sum(yes)-sum(z),6) for z in product(*legal)}
assert gaps == {F(-1,6),F(1,6)}
assert np.isclose(cov(np.ones(3),np.array(rows)/4,np.array(yes)/4),0)
x=np.array([0.,1.])
assert np.isclose(within(x,x,x),0) and np.all(x*(1-x)==0)
print("M2/M3: integer gaps -1/6,+1/6; CR empty. Pure units: both branches.")

# M4. Two 2x2 units, complements are the second column. Opposite column
# within signs jointly force D=ER. At T=.5 both individual cuts are feasible.
x = np.array([.25,.75])
# m0=(.5-x*m1)/(1-x); at these designs D = weights@m1 - 1.
weights = 2*x
def ends(A,b,bounds,c):
    opt = lambda s: linprog(s*c,A_ub=A,b_ub=b,bounds=bounds,method="highs")
    lo,hi=opt(1),opt(-1)
    assert lo.success and hi.success
    return np.array([lo.fun,-hi.fun])
bounds=[(0,1),(1/3,2/3)]
one=ends([weights],[1],bounds,weights)-1
both=ends([weights,-weights],[1,-1],bounds,weights)-1
assert np.allclose(one,[-.5,0]) and np.allclose(both,[0,0])
print("M4: joint within signs change [-.5,0] to {0}")

# M5. Two district intercepts and positive slopes in both rows.
x=np.array([.1,.2,.3,.7,.8,.9])
m_base=np.where(x<.5,.8,.2)+.1*x
aff=lambda k: within(x,m_base+k*(1-x),m_base-k*x)
assert aff(0)<0 and np.isclose(aff(.02),(aff(0)+aff(.04))/2)
print("M5: district-pooled within =",aff(0),"; affine in k")

# N2.1/N2.2. Extreme totals force both rates to be the total.
b=np.array([[0,0],[1,1]])
assert np.allclose(b.mean(0),[.5,.5]) and np.isclose(np.mean([.75,.25]),.5)
x=np.array([.25,.5,.75])
tau=x*.4+(1-x)*np.array([8/15,.7,.8])
assert np.allclose(tau,[.5,.55,.5])
print("N2: middle composition prevents constant endpoint extension")

# N2.3. A polynomial vanishing on the design changes the derivative.
z=Polynomial([0,1])
h=(z-.25)*(z-.5)*(z-.75)
assert np.allclose(h(x),0) and np.isclose(h.deriv()(.5),-1/16)
# Uniform X on [1/5,2/5]. Integrate each size-weighted covariance.
def E(f):
    integ=f.integ()
    return 5*(integ(.4)-integ(.2))
def weighted_term(w,m):
    return E(w*z*m)-E(w*z)*E(w*m)/E(w)
m1=z-z*z/2
m0=1-z*z/2
assert np.allclose((z*m1.deriv()+(1-z)*m0.deriv()).coef,0)
delta=weighted_term(z,m1)+weighted_term(1-z,m0)
assert np.isclose(delta,-1/47250,atol=1e-14)
print("N2.3: same support values, derivative shift -1/16; smooth local zero, aggregate",delta)
# Section 6. Check row-specific bin averages.
x=np.array([.25,.75])
assert np.isclose(np.mean(x*x),.3125) and np.isclose(x.mean()**2,.25)
assert np.isclose((x@x)/x.sum(),.625)

# R3(a,b). Full-law T in {0,1} fixes both means; T'=1/2 allows skew k.
x=np.array([.25,.5,.75])
for k in [.2,.3,.4]:
    m1,m0=.5+k*(1-x),.5-k*x
    assert np.allclose(x*m1+(1-x)*m0,.5)
    assert np.isclose(x@m1/x.sum()-(1-x)@m0/(1-x).sum(),5*k/6)
# R3(c). Two separated vertical strips, CR pair (u+.5,u-.5).
for target, expected in [(np.array([1,0]),[[-1,-.5],[.5,1]]),
                         (np.array([0,1]),[[-1,1],[-1,1]])]:
    pieces=[ends([[s,0]],[-.5],[(-1,1)]*2,target) for s in [1,-1]]
    assert np.allclose(pieces,expected)
# R3(d,e). Forced skew because off-diagonal probabilities sum to zero.
x=np.array([0.,.25,1.])
m1,m0=(49/55)*x,1-x
tau=x*m1+(1-x)*m0
d=x@m1/x.sum()-(1-x)@m0/(1-x).sum()
assert abs(cov(np.ones(3),x,tau))<1e-14 and np.isclose(d,-261/1925)
print("R3: sharp empty vs outer [1/6,1/3]; CR gap depends on target; LCE excludes CR")

# R2. Optimal facet with target intercept 2 and cuts that keep/remove it.
x_target=np.array([1,0])
for rhs,expect in [(1.2,3),(.8,2.8)]:
    assert np.isclose(2+ends([[1,1]],[rhs],[(0,1)]*2,x_target)[1],expect)
# Zero last-added gain for x, positive last-added gain for x+y.
for c in [x_target,np.ones(2)]:
    u0=ends([[1,0]],[.5],[(0,1)]*2,c)[1]
    u1=ends([[1,0],[1,1]],[.5,.8],[(0,1)]*2,c)[1]
    assert np.isclose(u0-u1,0 if c[1]==0 else .7)
print("R2: partial facet survives; full facet removed; target-specific redundancy")

# R10(a-d). Simulate an exact Gaussian experiment at ties, including emptiness.
rng=np.random.default_rng(831)
z1,z2,z3=rng.standard_normal((3,200000))
from scipy.optimize import brentq
q=brentq(lambda u:norm.cdf(u)**2*(2*norm.cdf(u)-1)-.95,1,4)
# k >= beta1 and k <= beta2; beta1=beta2=0, target k+beta3.
lo,hi=z1-q,z2+q
Egood=(z1<=q)&(-z2<=q)&(np.abs(z3)<=q)
cover=(lo<=hi)&(lo+z3-q<=0)&(0<=hi+z3+q)
assert np.all(cover[Egood])
# Parallel row k <= max(beta1,beta2); add k>=0, true k=0.
slab=np.maximum(z1+q,z2+q)
assert np.all(slab[(z1>=-q)&(z2>=-q)]>=0)
# W1 toy slab at beta=0; plus and minus branches calibrated separately
# with +/-Z directions, critical value z_.975. Empty if the ends cross.
t=norm.ppf(.975)
gamma=.25
pl,pu=gamma*(z1-t),z1+t
ml,mu=z1-t,gamma*(z1+t)
union=((pl<=0)&(0<=pu))|((ml<=0)&(0<=mu))
assert np.array_equal(union,abs(z1)<=t)
# Away from the tie, scaled true beta=4 and k=4 is on the true plus edge.
v=z1+4
cp=(gamma*(v-t)<=4)&(4<=v+t)
cm=(v-t<=4)&(4<=gamma*(v+t))
assert np.all((cp|cm)[cp])
print("R10: simultaneous, cover, empty, tie union, positive-edge union:",
      *[round(float(v),5) for v in [Egood.mean(),cover.mean(),(lo>hi).mean(),
                                     union.mean(),(cp|cm).mean()]])
# R10(e). Symmetric fixed design, fixed population edge, eight base rows
# and two target tails. The marginal-quantile union bound is deterministic.
x=np.array([.25,.5,.75])
design=np.array([(1-x)**2,2*x*(1-x),x*x]).T
p=np.ones(3)/3
er=np.array([cov(p,x,f)/cov(p,x,x) for f in design.T])
g=np.array([-((1-x)**2).sum()/(1-x).sum(),0,(x*x).sum()/x.sum()])
assert np.allclose(er,[-1,0,1]) and np.allclose(g,(7/12)*er)
print("R10(e): limiting naive coverage at most",norm.cdf((7/12)*norm.ppf(.995)))
from scipy.stats import binom
finite_draw_upper = norm.cdf((7/12)*norm.ppf(.996)) + binom.sf(199,4000,.04)
assert finite_draw_upper < .95
print("R10(e): bound with 4000 critical-value draws",finite_draw_upper)
