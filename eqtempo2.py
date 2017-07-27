# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 22:17:56 2016

@author: bapossatto
"""

import math
import numpy as np
import scipy
import matplotlib as mpl
from pylab import *

latitude=-28
longitude=-55.5
longitude_ref=-45.0

N=range(365) # Número de dias do ano
B=np.zeros(365) # vetor auxiliar
Gamma=np.zeros(365) # vetor auxiliar
E=np.zeros(365) # vetor armazenador dos resultados da equação do tempo

t_loc= np.repeat(np.array([range(24)]),365,axis=0) # matriz contendo o tempo local anual
t_sol=np.zeros(8760).reshape((365,24)) # matriz para os resultados do tempo solar aparente
omega=np.zeros(8760).reshape((365,24)) # matriz para os resultados do angulo horário
delta=np.zeros(8760).reshape((365,24)) # Ângulo de declinação solar
alpha=np.zeros(8760).reshape((365,24)) # Altitude solar
PSI=np.zeros(8760).reshape((365,24)) # Ângulo de azimute
PSI_n=np.zeros(8760).reshape((365,24)) # Ângulo de azimute
theta=np.zeros(8760).reshape((365,24)) #Ângulo de incidência
theta_z=np.zeros(8760).reshape((365,24)) #Ângulo complementar ao azimute
gamma=np.zeros(8760).reshape((365,24))
beta=np.zeros(8760).reshape((365,24))

for i in N:
 B[i] = (2*math.pi*(N[i]))/(365.0)
# B[i]=(N[i]-80)*360/364 
# ET[i]=9.87*math.sin(2*B[i]) - 7.53*math.cos(B[i]) - 1.5*math.sin(B[i])
 E[i] = 229.18*(0.000075 + 0.001868*math.cos(B[i]) - 0.032077*math.sin(B[i]) - 0.014615*math.cos(2*B[i]) - 0.04089*math.sin(2*B[i]))
 for ii in range(24):
     t_sol[i,ii] = t_loc[i,ii] + ((longitude_ref-longitude)/15.0) + E[i]/60.0
     omega[i,ii]=(12.0-t_sol[i,ii])*15*math.pi/180
     if sign(latitude)==1:
         delta[i,ii]= (0.006918  - 0.399912*math.cos(B[i]) + 0.070257*math.sin(B[i]) - 0.006758*math.cos(2*B[i]) + 0.000907*math.sin(2*B[i])  - 0.002697*math.cos(3*B[i]) + 0.00148*math.sin(3*B[i]))
     else:
         delta[i,ii]= (0.006918  - 0.399912*math.cos(B[i]) + 0.070257*math.sin(B[i]) - 0.006758*math.cos(2*B[i]) + 0.000907*math.sin(2*B[i])  - 0.002697*math.cos(3*B[i]) + 0.00148*math.sin(3*B[i]))
     alpha[i,ii] = math.asin(math.sin(latitude*math.pi/180)*math.sin(delta[i,ii]) + math.cos(latitude*math.pi/180)*math.cos(delta[i,ii])*math.cos(omega[i,ii]))
     if sign(omega[i,ii]) == 1:
         PSI[i,ii]=math.acos(((math.sin(alpha[i,ii])*math.sin(latitude*math.pi/180) - math.sin(delta[i,ii]))) / (math.cos(alpha[i,ii])*math.cos(latitude*math.pi/180)))
     else:
#         PSI[i,ii]=math.asin((math.cos(delta[i,ii])*math.sin(omega[i,ii]))/(math.cos(alpha[i,ii])))
         PSI[i,ii]= 2*math.pi - math.acos(((math.sin(alpha[i,ii])*math.sin(latitude*math.pi/180) - math.sin(delta[i,ii]))) / (math.cos(alpha[i,ii])*math.cos(latitude*math.pi/180)))
     theta[i,ii]=math.acos(math.sqrt(1-(math.cos(alpha[i,ii]-beta[i,ii])-math.cos(beta[i,ii])*math.cos(alpha[i,ii])*(1-math.cos(PSI[i,ii]-gamma[i,ii])))**2))
     PSI_n[i,ii]=(math.pi) - PSI[i,ii]
     theta_z[i,ii]=(math.pi/2)-alpha[i,ii]

xlabel('Tempo [dias]')
ylabel('Equacao do tempo [minutos]')
#title('$Performance$ $do$ $Veiulo$ $Eletrico$\n $Simulacao$ $da$ $aceleracao$')
plot(E, '+c')
#legend('$V_{max}$ $=$ $100$ $km/h$')
savefig('eq_tempo.png')
show()

xlabel('dia do ano')
ylabel('declinacao')
plot(N,delta, '-')
#plot(PSI*180/math.pi, -alpha*180/math.pi, '-b')
savefig('analema.png')
show()

xlabel('Psi_n')
ylabel('alpha')
plot(PSI_n*180/math.pi, alpha*180/math.pi, '-')
#plot(PSI*180/math.pi, -alpha*180/math.pi, '-b')
savefig('analema.png')
show()

alpha_novo=np.transpose(alpha)
PSI_NOVO=np.transpose(PSI)
plot(PSI_NOVO[17], alpha_novo[17],PSI_NOVO[16], alpha_novo[16],PSI_NOVO[10], alpha_novo[10],PSI_NOVO[9], alpha_novo[9],PSI_NOVO[5], alpha_novo[5])
