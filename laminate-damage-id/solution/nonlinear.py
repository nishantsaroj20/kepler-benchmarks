"""Condensed von Karman shallow-plate energy and implicit tangent sensitivities."""
import sys,json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,solve
from scipy.optimize import minimize,root
from numpy.polynomial.legendre import leggauss
from core import Linear,laminate

def extensional(s):
 # Separate thickness weighting; stiffness transformation is polynomial in ply angle.
 result=[]
 for group in [0,1]:
  a=np.zeros((3,3))
  for angle in s['layup_deg']:
   if int(abs(angle)==45)!=group:continue
   t=dict(s);t['layup_deg']=[angle];t['ply_thickness_m']=1.
   a+=sum(laminate(t))*12*s['ply_thickness_m']
  result.append(a)
 return np.array(result)

class Nonlinear:
 def __init__(self,s,pts,quad=(30,26)):
  self.s=s;self.L=s['Lx_m'];self.W=s['Ly_m'];self.h=s['ply_thickness_m']*len(s['layup_deg']);s=dict(s);s['model']=dict(s['model'],sine_orders=[3,3],gauss_orders=list(quad));self.lin=Linear(s,pts);self.n=9;self.mass=self.lin.mass
  self.xy=self.lin.xy;x,y=self.xy.T;al=self.lin.al;be=self.lin.be
  X=self.h*al*np.cos(x[:,None]*al)*np.sin(y[:,None]*be);Y=self.h*be*np.sin(x[:,None]*al)*np.cos(y[:,None]*be)
  self.imp=np.array(s['imperfection_h_coefficients']);gx=X@self.imp;gy=Y@self.imp
  nx,ny=quad;_,wx=leggauss(nx);_,wy=leggauss(ny);self.w=np.outer(wx,wy).ravel()*self.L*self.W/4
  self.B=np.zeros((len(x),3,18));self.B[:,0,:9]=X/self.h;self.B[:,1,9:]=Y/self.h;self.B[:,2,:9]=Y/self.h;self.B[:,2,9:]=X/self.h
  self.pairs=[(i,j) for i in range(9) for j in range(i,9)];nf=9+len(self.pairs);R=np.zeros((len(x),3,nf));R[:,0,:9]=gx[:,None]*X;R[:,1,:9]=gy[:,None]*Y;R[:,2,:9]=gx[:,None]*Y+gy[:,None]*X
  self.T2=np.zeros((nf,9,9))
  for k,(i,j) in enumerate(self.pairs,9):
   if i==j:R[:,0,k]=.5*X[:,i]**2;R[:,1,k]=.5*Y[:,i]**2;R[:,2,k]=X[:,i]*Y[:,i];self.T2[k,i,i]=2
   else:R[:,0,k]=X[:,i]*X[:,j];R[:,1,k]=Y[:,i]*Y[:,j];R[:,2,k]=X[:,i]*Y[:,j]+X[:,j]*Y[:,i];self.T2[k,i,j]=self.T2[k,j,i]=1
  self.R=R;self.A=extensional(s);self.omega=1000.;self.E0=self.mass*self.h**2*self.omega**2
  self.G=np.einsum('qi,qj,q->ij',X,X,self.w)+.35*np.einsum('qi,qj,q->ij',Y,Y,self.w)
  # Compression is specified once from the nominal undamaged model, then held fixed.
  self.N=float(s['precompression_N_m']);self.l=self.N*np.einsum('q,qi,q->i',gx,X,self.w)+.35*self.N*np.einsum('q,qi,q->i',gy,Y,self.w)
  self.bias=self.l/self.E0
  self.force=self.h*np.array([4/(m*n*np.pi**2) if m%2 and n%2 else 0 for m,n in self.lin.mn])/self.E0
 def assemble(self,p,deriv=False):
  gs=[self.lin.field(p[7*j:7*j+7]) for j in range((len(p)-2)//7)];f=np.prod([1-p[7*j+5:7*j+7,None]*gs[j] for j in range((len(p)-2)//7)],axis=0)
  a=np.einsum('gq,gab,q->qab',f,self.A,self.w)
  C=np.einsum('qai,qab,qbj->ij',self.B,a,self.B,optimize=True);D=np.einsum('qai,qab,qbj->ij',self.B,a,self.R,optimize=True)
  S=np.einsum('qai,qab,qbj->ij',self.R,a,self.R,optimize=True)-D.T@solve(C,D,assume_a='pos')
  self.S=S/self.E0;self.K=(self.lin.stiffness(p)*self.mass*self.h**2-self.N*self.G)/self.E0;self.bias=self.l/self.E0
  if deriv:
   gs2=[self.lin.field(p[7*j:7*j+7],True) for j in range((len(p)-2)//7)]
   factors=[1-p[7*j+5:7*j+7,None]*gs2[j][0] for j in range(len(gs2))];df=[]
   for j,(g,dg) in enumerate(gs2):
    d=p[7*j+5:7*j+7];other=np.prod([factors[k] for k in range(len(gs2)) if k!=j],axis=0)
    z=-d[:,None,None]*dg[None,:,:];z=np.concatenate([z,np.array([[-g,np.zeros_like(g)],[np.zeros_like(g),-g]]).transpose(0,2,1)],axis=2)
    df.extend((z*other[:,:,None]).transpose(2,0,1))
   df=np.array(df+[np.zeros_like(f)]*2);da=np.einsum('pgq,gab,q->pqab',df,self.A,self.w)
   Z=self.R-np.einsum('qai,ij->qaj',self.B,solve(C,D,assume_a='pos'))
   self.dS=np.einsum('qai,pqab,qbj->pij',Z,da,Z,optimize=True)/self.E0
   self.dK=self.lin.stiffness(p,True)[1]*self.mass*self.h**2/self.E0
  return self
 def damage_membrane(self, separated, loss):
  """Scale the intact membrane stiffness by 1-loss*s and rebuild the condensation."""
  At=self.A.sum(axis=0)
  factor=(1.0-float(loss)*np.asarray(separated,float))*self.w
  a=At[None,:,:]*factor[:,None,None]
  C=np.einsum('qai,qab,qbj->ij',self.B,a,self.B,optimize=True)
  D=np.einsum('qai,qab,qbj->ij',self.B,a,self.R,optimize=True)
  S=np.einsum('qai,qab,qbj->ij',self.R,a,self.R,optimize=True)-D.T@solve(C,D,assume_a='pos')
  self.S=S/self.E0
 def features(self,q):
  t=np.r_[q,[q[i]*q[j] for i,j in self.pairs]];dt=np.vstack([np.eye(9),[(np.eye(9)[i]*q[j]+np.eye(9)[j]*q[i]) for i,j in self.pairs]])
  return t,dt
 def energy(self,q,F):
  t,J=self.features(q);z=self.S@t;force=self.bias+F*self.force
  E=.5*q@self.K@q+.5*t@z-force@q;g=self.K@q+J.T@z-force;H=self.K+J.T@self.S@J+np.einsum('k,kij->ij',z,self.T2)
  return E,g,H
 def state(self,F,branch,start=None):
  q=np.zeros(9);q[0]=branch*1.5
  if start is not None:q=start
  sol=minimize(lambda v:self.energy(v,F)[:2],q,jac=True,hess=lambda v:self.energy(v,F)[2],method='trust-exact',options={'gtol':1e-10,'maxiter':100})
  q=sol.x
  for _ in range(5):
   E,g,H=self.energy(q,F)
   if np.linalg.norm(g)<1e-13:break
   if np.linalg.eigvalsh(H)[0]<=0:break
   q=q-solve(H,g,assume_a='pos')
  E,g,H=self.energy(q,F)
  if np.linalg.norm(g)>1e-7 or np.linalg.eigvalsh(H)[0]<=0 or np.sign(q[0])!=branch:raise ValueError('branch unavailable')
  return q,H
 def observe(self,F,branch,deriv=False):
  q,H=self.state(F,branch);val,V=eigh(H);fr=np.sqrt(val)*self.omega/(2*np.pi);Y=self.lin.H@V;Y/=np.linalg.norm(Y,axis=0)
  if not deriv:return q,fr,Y
  t,J=self.features(q);gp=np.einsum('pij,j->pi',self.dK,q)+np.einsum('ai,pab,b->pi',J,self.dS,t)
  dq=-solve(H,gp.T,assume_a='pos').T;dJ=np.einsum('kij,pj->pki',self.T2,dq)
  dH=self.dK+np.einsum('pai,ab,bj->pij',dJ,self.S,J)+np.einsum('ai,ab,pbj->pij',J,self.S,dJ)+np.einsum('ai,pab,bj->pij',J,self.dS,J)
  dH+=np.einsum('kij,pk->pij',self.T2,np.einsum('pab,b->pa',self.dS,t)+np.einsum('ab,bj,pj->pa',self.S,J,dq))
  P=np.einsum('ia,pij,jb->pab',V,dH,V,optimize=True);df=np.diagonal(P,axis1=1,axis2=2)*self.omega**2/(8*np.pi*np.pi*fr)
  gap=val[None,:]-val[:,None];np.fill_diagonal(gap,np.inf);dV=np.einsum('ia,pab->pib',V,P/gap)
  raw=self.lin.H@V;norm=np.linalg.norm(raw,axis=0);dY=np.einsum('si,pij->psj',self.lin.H,dV);dY=(dY-Y[None,:,:]*np.sum(Y[None,:,:]*dY,axis=1)[:,None,:])/norm[None,None,:]
  return q,fr,Y,dq,df,dY
 def fold(self,branch=-1):
  # Augmented equilibrium + null tangent vector, unknown transverse load.
  q,H=self.state(0,branch);val,V=eigh(H);v=V[:,0]
  def eq(z):
   q=z[:9];v=z[9:18];F=z[-1];_,g,H=self.energy(q,F);return np.r_[g,H@v,v@v-1]
  sol=root(eq,np.r_[q*.6,v,3.],tol=1e-10)
  if not sol.success or np.linalg.norm(eq(sol.x))>1e-7:raise ValueError(('fold failed',sol.message))
  return sol.x[-1],sol.x[:9],np.linalg.eigvalsh(self.energy(sol.x[:9],sol.x[-1])[2])

