# -*- coding: utf-8 -*-
"""
@author: bapossatto

No ciclo Rankine regenerativo mostrado na figura ao lado, uma vazão mássica de 89,42 kg/s de vapor a 16000 kPa e 540 o C é enviada da caldeira para a turbina. Após uma expansão parcial na turbina, uma fração de vapor y é extraída a 3200 kPa e enviada para um aquecedor de água de alimentação aberto, enquanto o restante do vapor segue expandindo na turbina até ser encaminhado a 45 o C ao condensador. Considerando que a eficiência isentrópica da turbina é 0,88 e que a eficiência isentrópica das bombas de alta e de baixa pressão é 0,85, pede-se:

a) Determinar a potência mecânica consumida pelas bombas, a potência térmica recebida na caldeira, a potência mecânica produzida pela turbina, a potência térmica rejeitada pelo condensador, a potência mecânica líquida do ciclo, a eficiência térmica do ciclo e a fração de vapor y enviada da turbina ao aquecedor.

b) Por que a eficiência térmica é menor que a calculada pelo arquivo “Ciclo Rankine ideal.EES”?

c) Obtenha um gráfico da eficiência térmica em função da eficiência isentrópica da turbina, fazendo-a variar de 0,7 a 1,0.

d) Obtenha um gráfico da eficiência térmica em função da eficiência isentrópica das bombas, fazendo-a variar de 0,7 a 1,0.

e) Obtenha um gráfico da eficiência térmica em função da pressão de extração do vapor na turbina, fazendo-a variar de 100 a 15000 kPa.}
"""

from iapws import IAPWS97
# from iapws import IAPWS95

E,E_s,P,h,T,v,x={},{},{},{},{},{},{} # [-],[-],[MPa],[kJ/kg],[K],[m^3/kg],[-]
s,s_r,h_r,Phase={},{},{},{} # [kJ/kg*K],[kJ/kg*K],[kJ/kg],[-]

m_dot=89.42 # [kg/s]
T[5]=540+273.15 # [K]
P[5]=16.0 # [MPa]
P[6]=3.2 # [MPa]
T[7]=45+273.15 # [K]
eta_s_t=0.88 # [-]
eta_s_b=0.85 # [-]

T[1]=T[7]
P[3]=P[6]
P[2]=P[3]
P[4]=P[5]

"1"
x[1]=0
E[1]=IAPWS97(T=T[1], x=x[1])
h[1],s[1],P[1],v[1],Phase[1]=E[1].h,E[1].s,E[1].P,E[1].v,E[1].region

"2"
s[2]=s[1]
E[2]=IAPWS97(P=P[2], s=s[2])
h[2],T[2],x[2],v[2],Phase[2]=E[2].h,E[2].T,E[2].x,E[2].v,E[2].region
h_r[2]=(h[2]-h[1]+(eta_s_b*h[1]))/eta_s_b

"3"
x[3]=0
E[3]=IAPWS97(P=P[3], x=x[3])
h[3],T[3],s[3],v[3],Phase[3]=E[3].h,E[3].T,E[3].s,E[3].v,E[3].region

"4"
s[4]=s[3]
E[4]=IAPWS97(P=P[4], s=s[4])
h[4],T[4],x[4],v[4],Phase[4]=E[4].h,E[4].T,E[4].x,E[4].v,E[4].region
h_r[4]=(h[4]-h[3]+(eta_s_b*h[3]))/eta_s_b

"5"
E[5]=IAPWS97(P=P[5], T=T[5])
h[5],s[5],x[5],v[5],Phase[5]=E[5].h,E[5].s,E[5].x,E[5].v,E[5].region

"6"
s[6]=s[5]
E[6]=IAPWS97(P=P[6], s=s[6])
h[6],T[6],x[6],v[6],Phase[6]=E[6].h,E[6].T,E[6].x,E[6].v,E[6].region
h_r[6]=h[5]-(eta_s_t*(h[5]-h[6]))

"7"
s[7]=s[6]
P[7]=P[1]
x[7]=(s[7]-0.6381)/7.5275
E[7]=IAPWS97(x=x[7], P=P[7])
h[7],T[7],x[7],v[7],Phase[7]=E[7].h,E[7].T,E[7].x,E[7].v,E[7].region
h_r[7]=h[6]-(eta_s_t*(h[6]-h[7]))

"Potência consumida pelas bombas"
y=(h[3]-h[2])/(h[6]-h[2])

W_dot_bbP=(1-y)*(h_r[2]-h[1])*m_dot

W_dot_baP=(h_r[4]-h[3])*m_dot

W_dot_b=W_dot_bbP+W_dot_baP
 
"Potência térmica recebida na caldeira"
W_dot_h=(h[5]-h[4])*m_dot

"Potência mecânica produzida pela turbina"
W_dot_t=m_dot*h[5]-(m_dot*((y*h[6])+((1-y)*h[7])))

"Potência térmica rejeitada no condensador"
Q_dot_l=(h[1]-h[7])*m_dot*(1-y)
W_dot_liq_abs=Q_dot_l*(-1)

"Potência mecânica líquida do ciclo"
W_dot_liq=W_dot_h-W_dot_liq_abs

"Eficiência Térmica do ciclo"
eta_term_cic=1 - (W_dot_liq_abs/W_dot_h)
eta = 1 - ((-1*Q_dot_l)/W_dot_h)