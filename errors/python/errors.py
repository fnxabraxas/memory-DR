# -*- coding: utf-8 -*-


def calcall1(eps,eps1,eps2,b,beta,N,nS1):
    import scipy.sparse.linalg as lin
    SD=np.zeros((nS1**2))
    payM,coop=calcPAYM(nS1,eps,eps1,eps2,b)
    #print(payM)
    fixM=pfixM(payM,nS1**2,beta,N)
    #print(fixM)
    val,v=lin.eigs(fixM,k=1,which='LR'); v=np.transpose(v.real); SD=v/np.sum(v)
    LC=np.dot(SD,coop)
    return SD,LC
    
    
def calcPAYM(nS1,eps,eps1,eps2,b):
    import scipy.sparse.linalg as lin
    payM=np.zeros((nS1**2,nS1**2))
    coop=np.zeros((nS1**2))
    st=np.zeros((nS1**2,nS1))
    c=1.
    Rp=b-c
    Pp=0.
    Sp=-c
    Tp=b
    #Rp=1.
    #Pp=0.
    #Sp=-0.5
    #Tp=1.5	
    for i in range (0,nS1**2):
        st[i,:]=dec2bin(np.array(list(map(int,(np.binary_repr(i, width=nS1))))),eps,eps1,eps2)
    for i in range (0,nS1**2):
        for j in range (0,nS1**2):
            M=createM(np.vstack((st[i],st[j])),nS1)
            val,v=lin.eigs(np.transpose(M),k=1,which='LR'); v=v.real; v=v/np.sum(v)
            payM[i,j]=v[0]*Rp+v[1]*Sp+v[2]*Tp+v[3]*Pp	
            #print(v)
            if(i==j): coop[i]=v[0]+v[1]
    return payM, coop


def dec2bin(p,eps,eps1,eps2):
    pf=np.zeros((4))
    p=(1.-eps)*p+eps*(1.-p)    
    pf[0]=(1.-eps2)*(1.-eps1)*p[0]+(1.-eps2)*eps1*p[2]+eps2*(1.-eps1)*p[1]+eps2*eps1*p[3]
    pf[1]=(1.-eps2)*(1.-eps1)*p[1]+(1.-eps2)*eps1*p[3]+eps2*(1.-eps1)*p[0]+eps2*eps1*p[2]
    pf[2]=(1.-eps2)*(1.-eps1)*p[2]+(1.-eps2)*eps1*p[0]+eps2*(1.-eps1)*p[3]+eps2*eps1*p[1]
    pf[3]=(1.-eps2)*(1.-eps1)*p[3]+(1.-eps2)*eps1*p[1]+eps2*(1.-eps1)*p[2]+eps2*eps1*p[0]
    return pf
      

def createM(p,nS1):
    M=np.zeros((nS1,nS1))
    p=np.transpose(p)

    M[0,0]=p[0,0]*p[0,1]
    M[0,1]=p[0,0]*(1.-p[0,1])
    M[0,2]=(1.-p[0,0])*p[0,1]
    M[0,3]=(1.-p[0,0])*(1.-p[0,1])
    
    M[1,0]=p[1,0]*p[2,1]
    M[1,1]=p[1,0]*(1.-p[2,1])
    M[1,2]=(1.-p[1,0])*p[2,1]
    M[1,3]=(1.-p[1,0])*(1.-p[2,1])
    
    M[2,0]=p[2,0]*p[1,1]
    M[2,1]=p[2,0]*(1.-p[1,1])
    M[2,2]=(1.-p[2,0])*p[1,1]
    M[2,3]=(1.-p[2,0])*(1.-p[1,1])
    
    M[3,0]=p[3,0]*p[3,1]
    M[3,1]=p[3,0]*(1.-p[3,1])
    M[3,2]=(1.-p[3,0])*p[3,1]
    M[3,3]=(1.-p[3,0])*(1.-p[3,1])

    return M


def pfix1(payII,payIJ,payJJ,payJI,betaN1,N):
    r=payII+payJJ-payIJ-payJI
    s=-payII+N*payIJ-(N-1.)*payJJ
    s2=-payJJ+N*payJI-(N-1.)*payII
    alp=np.exp(-betaN1*s)
    alp2=np.exp(-betaN1*s2)
    gam=-betaN1*r
    suma=0.
    suma2=0.
    term=1.
    term2=1.
    for k in range(1,N-1+1):
        term*=alp*np.exp(gam*k)
        term2*=alp2*np.exp(gam*k)
        suma+=term
        suma2+=term2
    pfix=1./(1.+suma)
    pfixI=1./(1.+suma2)
    return pfix,pfixI


def pfixM(payM,Ns,beta,N):
#	i becomes j (one individual of j take over the residents i)
    betaN1=beta/(N-1.)
    fixM=np.zeros((Ns,Ns))
    for i in range(0,Ns):
        for j in range (i+1,Ns):
            pfix,pfixI=pfix1(payM[i,i],payM[i,j],payM[j,j],payM[j,i],betaN1,N)
            fixM[i,j]=pfix
            fixM[j,i]=pfixI
    for i in range(0,Ns):
        fixM[i,i]=1.-np.sum(fixM[:,i])
    return fixM



if __name__ == "__main__":
    import numpy as np

    nS1=4
    beta=0.1
    N=1000
    
    eps=0.01
    eps1=0.
    eps2=0.
    b=3.

    SD,coop=calcall1(eps,eps1,eps2,b,beta,N,nS1)  
    print(SD)
    print(coop)
    
    #print(SD)
    #print(coop)


