"""
@author: kivancaykac
"""
#in this code we will use partial pivoting to solve kirchoff's law
# %% Libraries
import numpy as np
from cmath import polar, phase
import matplotlib.pyplot as plt
from Lab04_SolveLinear import PartialPivot
save = True  # Save the plots?

# %% Constants
R1, R3, R5 = 1e3, 1e3, 1e3  # kΩ
R2, R4, R6 = 2e3, 2e3, 2e3  # kΩ
C1 = 1e-6  # μF
C2 = .5e-6  # μF
xplus = 3.  # V
w = 1e3  # rad*s^-1

# %% First Part

def voltage(x,omega,t):
    return x*np.exp(1j*omega*t)

A = np.array([[(1/R1 + 1/R4 + 1j*w*C1), -1j*w*C1, 0.],
              [-1j*w*C1, (1/R2 + 1/R5 + 1j*w*(C1+C2)), -1j*w*C2],
              [0., -1j*w*C2, (1/R3 + 1/R6 + 1j*w*C2)]], complex)  # 3x3 matrix
v = np.array([xplus/R1, xplus/R2, xplus/R3], complex)  # 3x1

x = PartialPivot(A,v)
volt = voltage(x, w, 0.)
abs_volt = np.abs(volt)
phase_volt = np.angle(volt)*180/np.pi
print("Without Inductor")
print("|V_1|: {v1:.6f} V,  Phase_1 {p1:.6f} degrees\n|V_2|: {v2:.6f} V,  \
Phase_2 {p2:.6f} degrees\n|V_3|: {v3:.6f} V,  Phase_3 {p3:.6f} degrees".format(
v1=abs_volt[0],v2=abs_volt[1],v3=abs_volt[2],p1=phase_volt[0],p2=phase_volt[1],
p3=phase_volt[2]))

# First Part Plotting
t = np.linspace(0.,0.015,500)

volt1 = voltage(x[0],w,t)
phasevolt1 = np.angle(volt1)*180/np.pi
volt2 = voltage(x[1],w,t)
phasevolt2 = np.angle(volt2)*180/np.pi
volt3 = voltage(x[2],w,t)
phasevolt3 = np.angle(volt3)*180/np.pi

plt.style.use("seaborn-whitegrid")
plt.plot(t, voltage(xplus,w,t).real, label=r"$V_{+}$")
plt.plot(t, volt1.real, label=r"$V_1$")
plt.plot(t, volt2.real, label=r"$V_2$")
plt.plot(t, volt3.real, label=r"$V_3$")
plt.title(r"Voltage vs $Time$, without Inductor")
plt.ylabel(r"$Voltage$ $(V)$")
plt.xlabel(r"$Time$ $(sec)$")
plt.legend()
if(save): plt.savefig("Q1_without Inductor.png", dpi=1200)
plt.show()

# %% replace the resistor R6 with an inductor L
R6 = 1j*w*2  # =L
A = np.array([[(1/R1 + 1/R4 + 1j*w*C1), -1j*w*C1, 0.],
              [-1j*w*C1, (1/R2 + 1/R5 + 1j*w*(C1+C2)), -1j*w*C2],
              [0., -1j*w*C2, (1/R3 + 1/R6 + 1j*w*C2)]], complex)  # 3x3 matrix

x = PartialPivot(A,v)
volt = voltage(x, w, 0.)
abs_volt = np.abs(volt)
phase_volt = np.angle(volt)*180/np.pi
print("\n\nWith Inductor")
print("|V_1|: {v1:.6f} V,  Phase_1 {p1:.6f} degrees\n|V_2|: {v2:.6f} V,  \
Phase_2 {p2:.6f} degrees\n|V_3|: {v3:.6f} V,  Phase_3 {p3:.6f} degrees".format(
v1=abs_volt[0],v2=abs_volt[1],v3=abs_volt[2],p1=phase_volt[0],p2=phase_volt[1],
p3=phase_volt[2]))


volt1in = voltage(x[0],w,t)
phasevolt1in = np.angle(volt1in)*180/np.pi
volt2in = voltage(x[1],w,t)
phasevolt2in = np.angle(volt2in)*180/np.pi
volt3in = voltage(x[2],w,t)
phasevolt3in = np.angle(volt3in)*180/np.pi

plt.style.use("seaborn-whitegrid")
plt.plot(t, voltage(xplus,w,t).real, label=r"$V_{+}$")
plt.plot(t, volt1in.real, label=r"$V_1$")
plt.plot(t, volt2in.real, label=r"$V_2$")
plt.plot(t, volt3in.real, label=r"$V_3$")
plt.title(r"Voltage vs $Time$, with Inductor")
plt.ylabel(r"$Voltage$ $(V)$")
plt.xlabel(r"$Time$ $(sec)$")
plt.legend()
if(save): plt.savefig("Q1_with Inductor.png", dpi=1200)
plt.show()