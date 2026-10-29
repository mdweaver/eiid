# Illustration for ASSUMPTION_MENU.md section F (EI difference-in-differences).
# 2x2 per period: rows veteran (share x) / non-veteran; outcome "yes" on suffrage.
# Same counties, x and N in both periods; non-veteran contextual slope unchanged.
# Shows: the DiD does not depend on the unidentified slope (point identification),
# is exact without noise, and is centred on the truth with noise.
# Run: python3 -I illustrate_did.py
# 2x2: rows vet (share x) / non-vet; outcome "yes" for suffrage. Same counties, same x, same N in both periods.
import numpy as np
rng=np.random.default_rng(0); p=400
x=rng.uniform(.05,.6,p); N=rng.integers(500,5000,p).astype(float)
def run(bv0,bv1,bn0,bn1,noise=0.03):
    bv=np.clip(bv0+bv1*x+noise*rng.normal(size=p),0,1); bn=np.clip(bn0+bn1*x+noise*rng.normal(size=p),0,1)
    return x*bv+(1-x)*bn, bv, bn
# pre: future-vets slightly more pro; post: vets shift +0.25, non-vets +0.05; NON-VET SLOPE STABLE (-0.2) across periods
Tpre,bvp,bnp=run(0.30,0.10,0.25,-0.20); Tpost,bvq,bnq=run(0.55,0.10,0.30,-0.20)
wv=N*x; wn=N*(1-x)
truth = (np.sum(wv*bvq)/wv.sum()-np.sum(wv*bvp)/wv.sum()) - (np.sum(wn*bnq)/wn.sum()-np.sum(wn*bnp)/wn.sum())
def fit(T):
    return np.linalg.lstsq(np.c_[np.ones(p),x,x**2],T,rcond=None)[0]
def B_vet(th,w1): w0,c1,d1=th; return np.sum(wv*(w0+c1+d1*x-w1*(1-x)))/wv.sum()
def B_non(th,w1): w0,c1,d1=th; return np.sum(wn*(w0+w1*x))/wn.sum()
tp,tq=fit(Tpre),fit(Tpost)
# identified DiD under "non-vet contextual slope unchanged": w1_post = w1_pre = w (any value) -> DiD independent of w
for w in (-0.5,-0.2,0.3):
    did=(B_vet(tq,w)-B_vet(tp,w))-(B_non(tq,w)-B_non(tp,w))
    print(f"assumed common w1={w:+.1f}: DiD = {did:.4f}")
print(f"true DiD (finite population) = {truth:.4f}")
for lab,th,bv,bn in (("pre",tp,bvp,bnp),("post",tq,bvq,bnq)):
    print(lab,"theta",np.round(th,4),"B_vet est/true",round(B_vet(th,-0.2),4),round(np.sum(wv*bv)/wv.sum(),4),
          "B_non est/true",round(B_non(th,-0.2),4),round(np.sum(wn*bn)/wn.sum(),4))
# Monte Carlo with noise: is the identified DiD centred on the truth?
est=[];tru=[]
for s in range(300):
    rng=np.random.default_rng(1000+s)
    Tpre,bvp,bnp=run(0.30,0.10,0.25,-0.20); Tpost,bvq,bnq=run(0.55,0.10,0.30,-0.20)
    tp,tq=fit(Tpre),fit(Tpost)
    est.append((B_vet(tq,0)-B_vet(tp,0))-(B_non(tq,0)-B_non(tp,0)))
    tru.append((np.sum(wv*bvq)-np.sum(wv*bvp))/wv.sum()-(np.sum(wn*bnq)-np.sum(wn*bnp))/wn.sum())
est,tru=np.array(est),np.array(tru)
print(f"MC: mean est {est.mean():.4f}, mean truth {tru.mean():.4f}, sd of est-truth {np.std(est-tru):.4f}")
