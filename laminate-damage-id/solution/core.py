"""Spectral weak-form Kirchhoff laminate model, with analytic eigensensitivities."""
import numpy as np
from numpy.polynomial.legendre import leggauss

KEYS = ('xc','yc','a','b','theta_deg','d_cross','d_angle')

def laminate(s):
    p=s['material']; E1,E2,G,nu=p['E1_Pa'],p['E2_Pa'],p['G12_Pa'],p['nu12']
    den=1-nu*nu*E2/E1; q11,q22,q12,q66=E1/den,E2/den,nu*E2/den,G
    ds=[np.zeros((3,3)),np.zeros((3,3))]; h=s['ply_thickness_m']*len(s['layup_deg'])
    for i,angle in enumerate(s['layup_deg']):
        c,t=np.cos(np.deg2rad(angle)),np.sin(np.deg2rad(angle)); c2,t2=c*c,t*t
        q=np.zeros((3,3))
        q[0,0]=q11*c**4+2*(q12+2*q66)*c2*t2+q22*t**4
        q[1,1]=q11*t**4+2*(q12+2*q66)*c2*t2+q22*c**4
        q[0,1]=q[1,0]=(q11+q22-4*q66)*c2*t2+q12*(c**4+t**4)
        q[0,2]=q[2,0]=(q11-q12-2*q66)*c**3*t+(q12-q22+2*q66)*c*t**3
        q[1,2]=q[2,1]=(q11-q12-2*q66)*c*t**3+(q12-q22+2*q66)*c**3*t
        q[2,2]=(q11+q22-2*q12-2*q66)*c2*t2+q66*(c**4+t**4)
        z0=-h/2+i*s['ply_thickness_m']; z1=z0+s['ply_thickness_m']
        ds[int(abs(angle)==45)]+=q*(z1**3-z0**3)/3
    return ds

class Linear:
    def __init__(self,s,sensors,quadrature=None):
        self.s=s; self.L=s['Lx_m']; self.W=s['Ly_m']; self.nm=s['model']['sine_orders'][0]; self.nn=s['model']['sine_orders'][1]
        self.mn=np.array([(m,n) for m in range(1,self.nm+1) for n in range(1,self.nn+1)])
        self.al=self.mn[:,0]*np.pi/self.L; self.be=self.mn[:,1]*np.pi/self.W
        rho=s['material']['rho_kg_m3']; h=s['ply_thickness_m']*len(s['layup_deg'])
        self.mass=rho*h*self.L*self.W/4
        nx,ny=quadrature or s['model']['gauss_orders']; qx,wx=leggauss(nx); qy,wy=leggauss(ny)
        X,Y=np.meshgrid((qx+1)*self.L/2,(qy+1)*self.W/2,indexing='ij'); self.xy=np.array([X.ravel(),Y.ravel()]).T
        w=np.outer(wx,wy).ravel()*self.L*self.W/4/self.mass
        x,y=self.xy.T; sn=np.sin(x[:,None]*self.al)*np.sin(y[:,None]*self.be); cs=np.cos(x[:,None]*self.al)*np.cos(y[:,None]*self.be)
        B=np.stack([sn*self.al**2,sn*self.be**2,-2*cs*self.al*self.be],axis=1)
        self.T=np.array([np.einsum('qai,ab,qbj,q->qij',B,d,B,w,optimize=True) for d in laminate(s)])
        self.tf=self.T.reshape(2,len(x),-1); self.K0=self.T.sum((0,1)); self.edge=[]
        # w is fixed to zero on all four edges; only normal slopes have springs.
        for axis,edge in [('x',0),('x',1),('y',0),('y',1)]:
            if axis=='x':
                sign=(-1.)**self.mn[:,0] if edge else np.ones(len(self.mn)); v=self.al*sign
                e=np.outer(v,v)*(self.mn[:,1,None]==self.mn[None,:,1])*self.W/2
            else:
                sign=(-1.)**self.mn[:,1] if edge else np.ones(len(self.mn)); v=self.be*sign
                e=np.outer(v,v)*(self.mn[:,0,None]==self.mn[None,:,0])*self.L/2
            key=axis+('L' if edge else '0'); self.edge.append(e*s['rotational_springs_N'][key]/self.mass)
        self.edge=np.array(self.edge)
        x,y=np.asarray(sensors).T
        self.H=np.sin(x[:,None]*self.al)*np.sin(y[:,None]*self.be)

    def field(self,z,deriv=False):
        xc,yc,a,b,theta,dc,da=z; c,t=np.cos(theta*np.pi/180),np.sin(theta*np.pi/180)
        dx,dy=self.xy[:,0]-xc,self.xy[:,1]-yc; u=c*dx+t*dy; v=-t*dx+c*dy
        g=np.exp(-(u/a)**2-(v/b)**2)
        if not deriv:return g
        dg=np.stack([g*(2*u*c/a**2-2*v*t/b**2),g*(2*u*t/a**2+2*v*c/b**2),g*2*u*u/a**3,g*2*v*v/b**3,g*(-2*u*v/a**2+2*u*v/b**2)*np.pi/180],axis=1)
        return g,dg

    def stiffness(self,p,deriv=False):
        if p is None: return self.K0+self.edge.sum(0)
        count=(len(p)-2)//7
        fields=[self.field(p[7*j:7*j+7],deriv) for j in range(count)]
        gs=[f[0] if deriv else f for f in fields]
        factors=[1-np.asarray(p[7*j+5:7*j+7])[:,None]*gs[j] for j in range(count)]
        total=np.prod(factors,axis=0)
        K=(total.reshape(2,1,-1)@self.tf).sum(0).reshape(self.K0.shape)
        K+=self.edge[0]*p[-2]+self.edge[1]+self.edge[2]+self.edge[3]*p[-1]
        if not deriv:return K
        dK=[]
        for j in range(count):
            g,dg=fields[j];d=np.asarray(p[7*j+5:7*j+7]);other=np.prod([factors[k] for k in range(count) if k!=j],axis=0)
            dr=-d[:,None,None]*dg[None,:,:]
            dr=np.concatenate([dr,np.array([[-g,np.zeros_like(g)],[np.zeros_like(g),-g]]).transpose(0,2,1)],axis=2)
            dK.extend(np.einsum('gqp,gq,gqn->pn',dr,other,self.tf,optimize=True).reshape(7,*K.shape))
        return K,np.array(dK+[self.edge[0],self.edge[3]])

