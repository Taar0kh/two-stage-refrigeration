!pip install CoolPropimport numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot

## Two-stage vapor-compression refrigeration machine
## Incomplete interstage cooling
## Two throttling valves in series
## Refrigeration heat rate (Q_cold) is 90 kW
## Temperature of evaporation TEV = -60 °C
## Temperature of condensation TCD = 30 °C
## Refrigerants areR134 and R23 instead of R717

working_fluid=['R134a' , 'R23']
T = np.zeros((17, 1))
p = np.zeros((17, 1))
s = np.zeros((17, 1))
h = np.zeros((17, 1))

 #High Cascade
T_CD_HC = 30 + 273.16 #temperature of high temperature evaporator
T[4] = T_CD_HC #[K]
T_EV_LC = -60 + 273.16 #temperature of low temperature evaporator
T_m = (T_EV_LC*T_CD_HC)**0.5 #The intercascade temperature [K]

DELTAT_CD_HC=5  #[K]
DELTAT_EV_HC=5  #[K]
DELTAT_CD_LC=5  #[K]
DELTAT_EV_LC=5  #[K]

T[1]=T_m - DELTAT_EV_HC  #The evaporation temperature of the high cascade
p[1]=CP.PropsSI('P', 'T', T[1], 'Q', 1, working_fluid[0]) #pressure(R134a,T=T[1],x=1)
s[1]=CP.PropsSI('S', 'T', T[1], 'Q', 1, working_fluid[0]) #entropy(R134a,P=p[1],x=1)
h[1]=CP.PropsSI('H', 'T', T[1], 'Q', 1, working_fluid[0]) #enthalpy(R134a,P=p[1],x=1)

T_EV_HC=T[1]  #[K]

T[2]=T[1] + DELTAT_CD_HC
p[2]=p[1]
s[2]=CP.PropsSI('S', 'T', T[2], 'P', p[2], working_fluid[0]) #entropy(R134a,P=p[2],T=T[2])
h[2]=CP.PropsSI('H', 'T', T[2], 'P', p[2], working_fluid[0]) #enthalpy(R134a,P=p[2],T=T[2])
nu_2=CP.PropsSI('D', 'T', T[2], 'P', p[2], working_fluid[0]) #volume(R134a,P=p[2],T=T[2])

T_1_HC=T[2]  #[K]

p[3]=CP.PropsSI('P', 'T', T[4], 'Q', 1,    working_fluid[0]) #pressure(R134a,T=T[4],x=1)
T[3]=CP.PropsSI('T', 'S', s[2], 'P', p[3], working_fluid[0]) #temperature(R134a,s=s[3],P=p[3] )
s[3]=s[2]
h[3]=CP.PropsSI('H', 'T', T[3], 'P', p[3], working_fluid[0]) #enthalpy(R134a,P=p[3],T=T[3])

T[4]=T_CD_HC
p[4]=p[3]
s[4]=CP.PropsSI('S', 'T', T[4], 'Q', 1, working_fluid[0]) #entropy(R134a,P=p[4],x=1)
h[4]=CP.PropsSI('H', 'T', T[4], 'Q', 1, working_fluid[0]) #enthalpy(R134a,T=T[4],x=1)
T_0=T[4]-DELTAT_CD_HC

T[5]=T[4]
p[5]=p[4]

s[5]=CP.PropsSI('S', 'T', T[5], 'Q', 0, working_fluid[0]) #entropy(R134a,P=p[5],x=0)
h[5]=CP.PropsSI('H', 'T', T[5], 'Q', 0, working_fluid[0]) #enthalpy(R134a,P=p[5],x=0)

T[6]=T[5] - DELTAT_EV_HC
p[6]=p[5]
s[6]=CP.PropsSI('S', 'T', T[6], 'Q', 0, working_fluid[0]) #entropy(R134a,P=p[6],T=T[6])
h[6]=CP.PropsSI('H', 'T', T[6], 'Q', 0, working_fluid[0]) #enthalpy(R134a,P=p[6],s=s[6])

T[7]=T[1]
p[7]=p[1]
s_1=CP.PropsSI('S', 'T', T[7], 'Q', 0, working_fluid[0]) #entropy(R134a,P=p[5],x=0)
s_2=CP.PropsSI('S', 'T', T[7], 'Q', 1, working_fluid[0]) #entropy(R134a,P=p[5],x=0)
x_7=(s[6]-s_1)/(s_2-s_1)
s[7]=CP.PropsSI('S', 'T', T[6], 'Q',x_7, working_fluid[0])
h[7]=h[6]

T[8]=T[7]
p[8]=p[1]
h[8]=h[1]
s[8]=s[1] #Low Cascade
T[9]=-60+273.16#[K]
p[9]=CP.PropsSI('P', 'T', T[9], 'Q', 1, working_fluid[1]) #pressure(R23,T=T[9],x=1)
s[9]=CP.PropsSI('S', 'T', T[9], 'Q', 1, working_fluid[1]) #entropy(R23,P=p[9],x=1)
h[9]=CP.PropsSI('H', 'P', p[9], 'Q', 1, working_fluid[1]) #enthalpy(R23,P=p[9],x=1)

T_EV_LC=T[9]
T_cold=T[9]+DELTAT_EV_LC

T[10]=T[9] + DELTAT_CD_LC
p[10]=p[9]
s[10]=CP.PropsSI('S', 'T', T[10], 'P', p[10], working_fluid[1]) #entropy(R23,P=p[10],T=T[10])
h[10]=CP.PropsSI('H', 'P', p[10], 'T', T[10], working_fluid[1]) #enthalpy(R23,P=p[10],T=T[10])
nu_12=CP.PropsSI('D', 'T', T[10], 'P', p[10], working_fluid[1]) #volume(R23,P=p[10],T=T[10])

ETLC=T_m + DELTAT_CD_LC #The evaporation temperature of the low cascade

p[11]=CP.PropsSI('P', 'T', ETLC, 'Q', 1, working_fluid[1]) #pressure(R23,T=T[12],x=1)
T[11]=CP.PropsSI('T', 'S', s[10], 'P', p[11], working_fluid[1]) #temperature(R23,s=s[11],P=p[11] )
s[11]=s[10]
h[11]=CP.PropsSI('H', 'P', p[11], 'T', T[11], working_fluid[1]) #enthalpy(R23,P=p[11],T=T[11])

T[12]=ETLC
p[12]=p[11]
s[12]=CP.PropsSI('S', 'P', p[12], 'Q', 1, working_fluid[1]) #entropy(R23,P=p[12],x=1)
h[12]=CP.PropsSI('H', 'P', p[12], 'Q', 1, working_fluid[1]) #enthalpy(R23,T=T[12],x=1)

T_CD_LC=T[12]+273.16  #[K]

T[13]=T[12]
p[13]=p[11]
s[13]=CP.PropsSI('S', 'P', p[13], 'Q', 0, working_fluid[1]) #entropy(R23,P=p[13],x=0)
h[13]=CP.PropsSI('H', 'P', p[13], 'Q', 0, working_fluid[1]) #enthalpy(R23,P=p[13],x=0)

T[14]=T[13] - DELTAT_EV_LC
p[14]=p[11]
s[14]=CP.PropsSI('S', 'P', p[14], 'T', T[14], working_fluid[1]) #entropy(R23,P=p[14],T=T[14])
h[14]=CP.PropsSI('H', 'P', p[14], 'S', s[14], working_fluid[1]) #enthalpy(R23,P=p[14],s=s[14])

T[15]=T[9]
p[15]=p[9]
s_11=CP.PropsSI('S', 'T', T[15], 'Q', 0, working_fluid[1]) #entropy(R23,P=p[5],x=0)
s_22=CP.PropsSI('S', 'T', T[15], 'Q', 1, working_fluid[1]) #entropy(R23,P=p[5],x=0)
x_17=(s[14]-s_1)/(s_2-s_1)
s[15]=CP.PropsSI('S', 'T', T[14], 'Q',x_7, working_fluid[1])
h[15]=h[14]

T[16]=T[9]
p[16]=p[9]
h[16]=h[9]
s[16]=s[9]########################Temperature Entropy Diagrams##########################
plt.figure(figsize=(20,7))
plt.subplot(121)
plt.plot(s[1:9]/1000,T[1:9]-273.16,"*-r")
plt.plot(s[9:17]/1000,T[9:17]-273.16,"*-b")
plt.xlabel('s [kJ/kg/K]')
plt.ylabel('T [C]')
plt.title("Temperature-Entropy Diagram")

##########################Bell Diagram##########################
Tcritp=CP.PropsSI("Tcrit",working_fluid[0])
Tmin_p=CP.PropsSI("T_MIN",working_fluid[0])
i= np.linspace(Tmin_p,Tcritp,1000)
sf=CP.PropsSI('S', 'T', i, 'Q', 0, working_fluid[0])
sv=CP.PropsSI('S', 'T', i, 'Q', 1, working_fluid[0])
plt.plot(sf/1000,i-273.16,'--r',linewidth=.5)
plt.plot(sv/1000,i-273.16,'--r',linewidth=.5)

# working_fluid=working_fluid[0]
tcritp=CP.PropsSI("Tcrit",working_fluid[1])
tmin_p=CP.PropsSI("T_MIN",working_fluid[1])
j= np.linspace(tmin_p,tcritp,1000)
sf=CP.PropsSI('S', 'T', j, 'Q', 0, working_fluid[1])
sv=CP.PropsSI('S', 'T', j, 'Q', 1, working_fluid[1])
plt.plot(sf/1000,j-273.16,'--b',linewidth=.5)
plt.plot(sv/1000,j-273.16,'--b',linewidth=.5)
################################################################

for i in range(1,17):
    if i != 8:
      plt.text(s[i]/1000,T[i]-273.16, i-8*(i>8),{'size':10})

plt.text(0.5,max(T)-263, '$R23$', {'color':'blue','size':20})
plt.text(0.5,max(T)-273, '$R134a$',{'color':'red','size':20})

plt.axis([0.5, 2, -75, 75])

########################Enthalpy Pressure Diagrams##########################
plt.subplot(122)
plt.plot(h[1:9]/1000,p[1:9]/1000,".-r")
plt.plot(h[9:17]/1000,p[9:17]/1000,".-b")
plt.yscale('log')
plt.xlabel('h [kJ/kg]')
plt.ylabel('p [kPa]')
plt.title("Pressure-Enthalpy Diagram")

##########################Bell diagram##########################
Pcritp=CP.PropsSI("Pcrit",working_fluid[0])
Pmin_p=CP.PropsSI("P_MIN",working_fluid[0])
j= np.linspace(Pmin_p,Pcritp,1000)
hf=CP.PropsSI('H', 'P', j, 'Q', 0, working_fluid[0])
hv=CP.PropsSI('H', 'P', j, 'Q', 1, working_fluid[0])
plt.plot(hf/1000,j/1000,'--r',linewidth=.5)
plt.plot(hv/1000,j/1000,'--r',linewidth=.5)


pcritp=CP.PropsSI("Pcrit",working_fluid[1])
pmin_p=CP.PropsSI("P_MIN",working_fluid[1])
j= np.linspace(pmin_p,pcritp,1000)
hf=CP.PropsSI('H', 'P', j, 'Q', 0, working_fluid[1])
hv=CP.PropsSI('H', 'P', j, 'Q', 1, working_fluid[1])
plt.plot(hf/1000,j/1000,'--b',linewidth=.5)
plt.plot(hv/1000,j/1000,'--b',linewidth=.5)
################################################################

for j in range(1,16):
  if j != 8:
    plt.text(h[j]/1000,p[j]/1000, j-8*(j>8),{'size':10})
plt.axis([150, 450, 50, 4000])
plt.show()#Data
c=0.03 #[]
m_1=1 #[]
m_2=1 #[]
eta_e_motor=0.9 #[]
eta_e_motor1=0.9 #[]
b=0.0025 #[]
p_i_fr_LC=40  #[kPa]
p_i_fr_HC=40  #[kPa]
alpha=1.12 #[]
beta=0.5 #[]
q_dot_cold_LC=90  #[kW]
T_1_LC=T[10]+273.16  #[K]
teta_LC=T_1_LC-T_EV_LC
q_cold_LC=h[10]-h[15]
q_v_LC=q_cold_LC/nu_12
q_CD_LC=h[11]-h[14]
m_dot_wf_LC=q_dot_cold_LC/q_cold_LC
q_dot_CD_LC=q_CD_LC*m_dot_wf_LC
q_dot_EV_HC=q_dot_CD_LC
teta_HC=T_1_HC-T_EV_HC
q_cold_HC=h[2]-h[7]
q_v_HC=q_cold_HC/nu_2
q_CD_HC=h[3]-h[6]
w_cyc_HC=h[3]-h[2]
m_dot_wf_HC=q_dot_EV_HC/q_cold_HC
V_dot_real_HC=m_dot_wf_HC*nu_2
lambda_c_HC=1-c*((p[3]/p[1])**(1/m_1)-1)
lambda_w_HC=(T_EV_HC+teta_HC)/(alpha*T_CD_HC+0.5*teta_HC)
lambda_HC=lambda_c_HC*lambda_w_HC
V_dot_h_HC=V_dot_real_HC/lambda_HC
eta_i_HC=lambda_w_HC+b*T[1]
w_cyc_LC=h[11]-h[10]
V_dot_real_LC=m_dot_wf_LC*nu_12
lambda_c_LC=1-c*((p[11]/p[9])**(1/m_2)-1)
lambda_w_LC=(T_EV_LC+teta_LC)/(alpha*T_CD_LC+0.5*teta_LC)
eta_i_LC=lambda_w_LC+b*T[9]
lambda_LC=lambda_c_LC*lambda_w_LC
V_dot_h_LC=V_dot_real_LC/lambda_LC
P_dot_th_LC=m_dot_wf_LC*w_cyc_LC
P_dot_th_HC=m_dot_wf_HC*w_cyc_HC
P_dot_i_LC=P_dot_th_LC/eta_i_LC
P_dot_i_HC=P_dot_th_HC/eta_i_HC
P_dot_fr_LC=V_dot_h_LC*p_i_fr_LC
P_dot_fr_HC=V_dot_h_HC*p_i_fr_HC
P_dot_e_LC=P_dot_i_LC+P_dot_fr_LC
P_dot_e_HC=P_dot_i_HC+P_dot_fr_HC
COP_star_carnot=(T_cold)/(T_0-T_cold)
COP_carnot=T_EV_LC/(T_CD_HC-T_EV_LC)

COP_th_RM=q_dot_cold_LC/(P_dot_th_LC+P_dot_th_HC)
COP_real_RM=q_dot_cold_LC/(P_dot_e_LC+P_dot_e_HC)

eta_star_RM=COP_th_RM/COP_star_carnot
eta_RM=COP_th_RM/COP_carnot
eta=COP_real_RM/COP_carnot

Name=['c', 'm_1', 'm_2', 'eta_e_motor', 'eta_e_motor1', 'b', 'p_i_fr_LC', 'p_i_fr_HC', 'alpha', 'beta', 'q_dot_cold_LC', 'T_1_LC', 'teta_LC', 'q_cold_LC', 'q_v_LC', 'q_CD_LC', 'm_dot_wf_LC', 'q_dot_CD_LC', 'q_dot_EV_HC', 'teta_HC', 'q_cold_HC', 'q_v_HC', 'q_CD_HC', 'w_cyc_HC', 'm_dot_wf_HC', 'V_dot_real_HC', 'lambda_c_HC', 'lambda_w_HC', 'lambda_HC', 'V_dot_h_HC', 'eta_i_HC', 'w_cyc_LC', 'V_dot_real_LC', 'lambda_c_LC', 'lambda_w_LC', 'eta_i_LC', 'lambda_LC', 'V_dot_h_LC', 'P_dot_th_LC', 'P_dot_th_HC', 'P_dot_i_LC', 'P_dot_i_HC', 'P_dot_fr_LC', 'P_dot_fr_HC', 'P_dot_e_LC', 'P_dot_e_HC', 'COP_star_carnot', 'COP_carnot', 'COP_th_RM', 'COP_real_RM', 'eta_star_RM', 'eta_RM', 'eta']
variables=[c, m_1, m_2, eta_e_motor, eta_e_motor1, b, p_i_fr_LC, p_i_fr_HC, alpha, beta, q_dot_cold_LC, T_1_LC, teta_LC, q_cold_LC, q_v_LC, q_CD_LC, m_dot_wf_LC, q_dot_CD_LC, q_dot_EV_HC, teta_HC, q_cold_HC, q_v_HC, q_CD_HC, w_cyc_HC, m_dot_wf_HC, V_dot_real_HC, lambda_c_HC, lambda_w_HC, lambda_HC, V_dot_h_HC, eta_i_HC, w_cyc_LC, V_dot_real_LC, lambda_c_LC, lambda_w_LC, eta_i_LC, lambda_LC, V_dot_h_LC, P_dot_th_LC, P_dot_th_HC, P_dot_i_LC, P_dot_i_HC, P_dot_fr_LC, P_dot_fr_HC, P_dot_e_LC, P_dot_e_HC, COP_star_carnot, COP_carnot, COP_th_RM, COP_real_RM, eta_star_RM, eta_RM, eta]
for i in range(1,len(Name)):
  print("%15s" % Name[i],' = ', float(np.round(variables[i],3))))
