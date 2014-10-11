#-*- coding: utf-8 -*-
import numpy as n
from scipy.io import wavfile as w

H=n.hstack
V=n.vstack

f_a = 44100. # Hz, frequência de amostragem

############## 2.2.1 Tabela de busca (LUT)
Lambda_tilde=Lt=1024.*16

# Senoide
foo=n.linspace(0,2*n.pi,Lt,endpoint=False)
S_i=n.sin(foo) # um período da senóide com T amostras

# Quadrada:
Q_i=n.hstack(  ( n.ones(Lt/2)*-1 , n.ones(Lt/2) )  )

# Triangular:
foo=n.linspace(-1,1,Lt/2,endpoint=False)
Tr_i=n.hstack(  ( foo , foo*-1 )   )

# Dente de Serra:
D_i=n.linspace(-1,1,Lt)

def v(f=200,d=2.,tab=S_i,fv=2.,nu=2.,tabv=S_i):
    Lambda=n.floor(f_a*d)
    ii=n.arange(Lambda)
    Lv=float(len(tabv))

    Gammav_i=n.floor(ii*fv*Lv/f_a) # índices para a LUT
    Gammav_i=n.array(Gammav_i,n.int)
    # padrão de variação do vibrato para cada amostra
    Tv_i=tabv[Gammav_i%int(Lv)] 

    # frequência em Hz em cada amostra
    F_i=f*(   2.**(  Tv_i*nu/12.  )   ) 
    # a movimentação na tabela por amostra
    D_gamma_i=F_i*(Lt/float(f_a))
    Gamma_i=n.cumsum(D_gamma_i) # a movimentação na tabela total
    Gamma_i=n.floor( Gamma_i) # já os índices
    Gamma_i=n.array( Gamma_i, dtype=n.int) # já os índices
    return tab[Gamma_i%int(Lt)] # busca dos índices na tabela

def A(fa=2.,V_dB=10.,d=2.,taba=S_i):
    Lambda=n.floor(f_a*d)
    ii=n.arange(Lambda)
    Lt=float(len(taba))
    Gammaa_i=n.floor(ii*fa*Lt/f_a) # índices para a LUT
    Gammaa_i=n.array(Gammaa_i,n.int)
    # variação da amplitude em cada amostra
    A_i=taba[Gammaa_i%int(Lt)] 
    A_i=A_i*10.**(V_dB/20.)
    return A_i



def adsr(som,A=10.,D=20.,S=-20.,R=100.,xi=1e-2):
    a_S=10**(S/20.)
    Lambda=len(som)
    Lambda_A=int(A*f_a*0.001)
    Lambda_D=int(D*f_a*0.001)
    Lambda_R=int(R*f_a*0.001)

    ii=n.arange(Lambda_A,dtype=n.float)
    A=ii/(Lambda_A-1)
    A_i=A
    ii=n.arange(Lambda_A,Lambda_D+Lambda_A,dtype=n.float)
    D=1-(1-a_S)*(   ( ii-Lambda_A )/( Lambda_D-1) )
    A_i=n.hstack(  (A_i, D  )   )
    S=n.ones(Lambda-Lambda_R-(Lambda_A+Lambda_D),dtype=n.float)*a_S
    A_i=n.hstack( ( A_i, S )  )
    ii=n.arange(Lambda-Lambda_R,Lambda,dtype=n.float)
    R=a_S-a_S*((ii-(Lambda-Lambda_R))/(Lambda_R-1))
    A_i=n.hstack(  (A_i,R)  )

    return som*A_i

triadeM=[0.,4.,7.]
def ac(f=200.,notas=[0.,4.,7.,12.],tab=Q_i):
    acorde=adsr(v(tab=tab,f=f*2.**(notas[-1]/12.),nu=0))
    for na in notas[:-1]:
        acorde+=adsr(v(tab=tab,f=f*2**(na/12.),nu=0))
    
    return acorde

s=ac(200.,triadeM)
s2=ac(200.,triadeM,Tr_i)
#s=n.hstack(( s,intervaloHarmonico(200,i) ))
s=H((s,s2,s,s2))


s1=ac(200.,triadeM)
s2=ac(200.,[0,5,9]) # subdominante
s3=ac(200.,[2,7,11]) # dominante
s4=ac(200.,[2,5,9]) # sub relativa
s5=ac(200.,[0,4,9]) # ton relat / sub anti


s=H((s,s1,s2,s3,s1,  s1,s4,s3,s1,  s1,s2,s3,s5, s5,s4,s3,s1  ))
s_=n.copy(s)


s=((s-s.min())/(s.max()-s.min()))*2.-1.

# most music players read only 16-bit wav files, so let's convert the array
s = n.int16(s * float(2**15-1))
w.write("acordeCedo.wav",f_a, s) # escrita do som

f0=100.
padrao=[0,0,-1,0]
padrao2=[10,0,-1,50]
padrao3=[10,0,-1,50]
padrao4=[10,0,-5,5]
padrao5=[11,0,-1,7]
pf=[f0*2**(pi/12.) for pi in padrao]
som1=n.hstack([adsr(v(f=ff,d=0.5,tab=Q_i)) for ff in pf])
som2=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao2])
som3=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao3])
som4=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao4])
som5=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao5])

s=n.hstack((som1,som2,som5,som1,som4,som2,som3,som1))
s__=n.copy(s)
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("linha.wav",f_a, s) # escrita do som


# para sobrepor com s_
# Intro:
padrao=[0,0,-1,0]
padrao1=[0,4,7,12]
padrao2=[0,-5,-9,-12]
padrao3=[-9,-5,-3,-1]
pp=padrao+padrao1+padrao2+padrao3
som=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in pp])

#Ciclo1
padrao1=[(i,.5) for i in [0,12,24,36]]
padrao2=[(0,.25),(0,.25),(0,.5),(24+5,.5),(36+9,.5)]
padrao3=[(0,.25),(-1,.25),(-5,.5),(24+7,.5),(36+11,.5)]
padrao4_1=[(0,.25),(0,.25),(0,.5),(-24,.5),(36*2,.5)]
pp_=padrao1+padrao2+padrao3+padrao4_1

#Ciclo2
padrao1=[(i,.5) for i in [0,12,24,36]]
padrao2=[(0,.25),(2,.25),(2,.5),(24+2,.5),(36+9,.5)]
padrao3=[(0,.25),(-1,.25),(-5,.5),(24-1,.5),(36+2,.5)]
padrao4_1=[(0,.25),(-5,.25),(-9,.5),(-12,.5),(-5,.5)]
pp_+=padrao1+padrao2+padrao3+padrao4_1

#Ciclo3
padrao1=[(i,.5) for i in [0,12,24,36]]
padrao2=[(0,.25),(0,.25),(0,.5),(24+2,.5),(36+9,.25),(36+12,.25)]
padrao3=[(11,.25),(-1,.25),(7,.5),(24+7,.5),(36+2,.25),(36-1,.25)]
padrao5=[(0,.25),(-3,.25),(4,.5),(9,.5),(-5,.5)]
padrao5b=[(0,.25),(-3,.25),(4,.25),(9,.25),(-5,.25),(0,.25),(-3,.25),(4,.25)]
padrao4=[(0,.25),(2,.25),(2,.5),(24+2,.5),(36+9,.25),(-36+9,.25)]
padrao3_=[(11,.25),(-1,.25),(7,.25),(24+7,.5),(36+2,.5),(36-1,.25)]
padrao=   [(i,.5) for i in [0,0,-1,0]]
padrao1_= [(i,.5) for i in [0,4,7,12]]
padrao2_= [(i,.5) for i in [0,-5,-9,-12]]
padrao3__=[(i,.5) for i in [-9,-5,-3,-1]]
padraoF=  [(i,.5) for i in [-12]]
pp_+=padrao1+padrao2+padrao3+padrao5+padrao5b+padrao4+padrao3_+padrao+padrao1_+padrao2_+padrao3__+padraoF

som_=n.hstack([adsr(v(f=f0*2**(pi[0]/12.),d=pi[1],tab=Q_i)) for pi in pp_])

som=n.hstack((som,som_))


# ve a diferenca enre s_ e som, adiciona zeros em s_ e soma ambos.
# coloca s__ no comeco de tudo. Grava wav
s=n.hstack((s__,   n.hstack((  s_,n.zeros(len(som)-len(s_))   ))+som ))
###s= n.hstack((  s_,n.zeros(len(som)-len(s_))   ))+som 

s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("lunhani.wav",f_a, s) # escrita do som
