import numpy as n
import matplotlib.pyplot as p
import scipy as s

# V_PxG  W_PxG

V000 = 8718 # P = 14, x = -0.5, G = 4500
V001 = 8978 # P = 14, x = -0.5, G = 5000
V010 = 7543 # P = 14, x = -0.4, G = 4500
V011 = 7649 # P = 14, x = -0.4, G = 5000
V100 = 8015 # P = 16, x = -0.5, G = 4500
V101 = 8249 # P = 16, x = -0.5, G = 5000
V110 = 6902 # P = 16, x = -0.4, G = 4500
V111 = 7179 # P = 16, x = -0.4, G = 5000

W000 = 7084 # P = 14, x = -0.5, G = 4000
W001 = 5339 # P = 14, x = -0.5, G = 5000
W010 = 6131 # P = 14, x = -0.4, G = 4000
W011 = 4781 # P = 14, x = -0.4, G = 5000
W100 = 5835 # P = 16, x = -0.5, G = 4000
W101 = 4467 # P = 16, x = -0.5, G = 5000
W110 = 4900 # P = 16, x = -0.4, G = 4000
W111 = 3805 # P = 16, x = -0.4, G = 5000

G0 = 4000
G1 = 4500
G2 = 5000

x1 = -0.5
x2 = -0.4

P1 = 14
P2 = 16

s = 12.75 # mm
d = 9.1 # mm
D_h = 10.5 # mm

x = [-0.428252, -0.427085, -0.425905, -0.424745, -0.423599, -0.422459, -0.42132, -0.420184, -0.419079, -0.418145]
G = [4513.528727783, 4504.715675237, 4498.613443827, 4494.3914841, 4491.473025599, 4489.419722781, 4488.081615326, 4486.881932781, 4486.974216054, 4483.859655599]
P = [15.5737, 15.5653, 15.5569, 15.5485, 15.54, 15.5316, 15.5232, 15.5148, 15.5063, 15.4979]

H = 3.53

#-------------------------------------------------------------------------------
# A Trilinear Interpolation Algorithm by Paul Bourke
# Link: http://paulbourke.net/miscellaneous/interpolation/

# a = scaled pressure
# b = scaled equilibrium quality
# c = scaled mass flux
# V_xyz = V_abc = V_PXG

# Scaling pressure

a = []
for i in range(10):
    z = (P[i] - 14)/2
    a.append(z)

# Scaling equilibrium quality

b = []
for i in range(10):
    z = (x[i] + 0.5)/0.1
    b.append(z)

# Scaling mass flux

c = []
for i in range(10):
    z = (G[i] - 4500)/500
    c.append(z)

CHF_G = []
for i in range(10):
    z = V000*(1 - a[i])*(1 - b[i])*(1 - c[i]) + V100*a[i]*(1 - b[i])*(1 - c[i]) + V010*(1 - a[i])*b[i]*(1 - c[i]) + V001*(1 - a[i])*(1 - b[i])*c[i] + V101*a[i]*(1 - b[i])*c[i] + V011*(1 - a[i])*b[i]*c[i] + V110*a[i]*b[i]*(1 - c[i]) + V111*a[i]*b[i]*c[i]
    CHF_G.append(z)

f = []
for i in range(10):
    z = (G[i] - 4000)/1000
    f.append(z)

CHF_B = []
for i in range(10):
    z = W000*(1 - a[i])*(1 - b[i])*(1 - f[i]) + W100*a[i]*(1 - b[i])*(1 - f[i]) + W010*(1 - a[i])*b[i]*(1 - f[i]) + W001*(1 - a[i])*(1 - b[i])*f[i] + W101*a[i]*(1 - b[i])*f[i] + W011*(1 - a[i])*b[i]*f[i] + W110*a[i]*b[i]*(1 - f[i]) + W111*a[i]*b[i]*f[i]
    CHF_B.append(z)

#-------------------------------------------------------------------------------
# Correction factors for Groeneveld look-up table

# Wong exponent
exponent = []
for i in range(10):
    a = 0.58*(1-0.25*n.exp(-2*x[i]))*(1-15*D_h**(-6)*G[i])
    exponent.append(a)

print(n.round(exponent, 4))
print(type(n.round(exponent, 4)))

K1_G1 = []
for i in range(10):
    a = (8/D_h)**(exponent[i])
    K1_G1.append(a)

# Groeneveld exponent
K1_G2 = []
for i in range(10):
    a = (8/D_h)**(0.5)
    K1_G2.append(a)

# Tanase exponent
K1_G3 = []
for i in range(10):
    a = (8/D_h)**(0.3)
    K1_G3.append(a)

K2_G = []
for i in range(10):
    a = (0.5+2*(s- d)/d)*n.exp(-0.5*(n.abs(x[i]))**(1/3))
    K2_G.append(a)

alpha = [1.792456951, 1.794026504, 1.796968338, 1.799242769, 1.801253564, 1.803205343, 1.805172303, 1.807145907, 1.809112296, 1.811060236]

K4_G = []
for i in range(10):
    L = 0.353*(i+1)
    j = n.exp((0.0105/L)*n.exp(2*alpha[i]))
    K4_G.append(j)

CHF_GI1 = [] # Wong
for i in range(10):
    g = K1_G1[i]*K2_G[i]*K4_G[i]*CHF_G[i]
    CHF_GI1.append(g)

CHF_GI2 = [] # Groeneveld
for i in range(10):
    g = K1_G2[i]*K2_G[i]*K4_G[i]*CHF_G[i]
    CHF_GI2.append(g)

CHF_GI3 = [] # Tanase
for i in range(10):
    g = K1_G3[i]*K2_G[i]*K4_G[i]*CHF_G[i]
    CHF_GI3.append(g)
#-------------------------------------------------------------------------------
# Correction factors for Bobkov look-up table

K1_B = (9.36/D_h)**(1/3)

K2_B = 0.2 + 0.57*s/d

K3_B = []
for i in range(10):
    L = 0.353*(i+1)
    a = 1 + 0.6 * n.exp(-0.01*L/0.0105)
    K3_B.append(a)

K4_B = [1+1.6*(G[i]/1000)*n.exp(-(H/15)/(7.6*D_h/1000))*(0.2+n.abs(x[i])) for i in range(10)]
print("F_B", K4_B)

CHF_BI = []
for i in range(10):
    a = K1_B*K2_B*K3_B[i]*K4_B[i]*CHF_B[i]
    CHF_BI.append(a)
#-------------------------------------------------------------------------------
# Error estimation for diameter factor exponent
# Normalized Mean Absolute Error (NMAE)
av = sum(CHF_BI)/10

Sum1 = 0
Sum2 = 0
Sum3 = 0
for i in range(10):
    Sum1 = Sum1 + n.abs(CHF_BI[i] - CHF_GI1[i])
    Sum2 = Sum2 + n.abs(CHF_BI[i] - CHF_GI2[i])
    Sum3 = Sum3 + n.abs(CHF_BI[i] - CHF_GI3[i])
ER1 = Sum1*100/10/av
ER2 = Sum2*100/10/av
ER3 = Sum3*100/10/av

# Normalized Root Mean Square Error (NRMSE)

Sum4 = 0
Sum5 = 0
Sum6 = 0
for i in range(10):
    Sum4 = Sum4 + n.abs(CHF_BI[i] - CHF_GI1[i])**2
    Sum5 = Sum5 + n.abs(CHF_BI[i] - CHF_GI2[i])**2
    Sum6 = Sum6 + n.abs(CHF_BI[i] - CHF_GI3[i])**2
ER4 = n.sqrt(Sum4/10)*100/av
ER5 = n.sqrt(Sum5/10)*100/av
ER6 = n.sqrt(Sum6/10)*100/av

CHFI   = n.array(CHF_GI1)
CHFII  = n.array(CHF_GI2)
CHFIII = n.array(CHF_GI3)
CHFIV  = n.array(CHF_BI)

N = n.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

p.plot(N, CHFI,   "--sb", label = "Wong exponent. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER1, ER4))
p.plot(N, CHFII,  "-.Dg", label = "Groeneveld exponent. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER2, ER5))
p.plot(N, CHFIII, "--*r", label = "Tanase exponent. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER3, ER6))
p.plot(N, CHFIV,  "--om", label = "Bobkov LUT")

p.grid("on")
p.xlabel("Normalized length ($z/H$)")
p.ylabel("CHF ($kW/m^{2}$)")
p.legend()
p.show()

# Empirical correlations

CHFW3 = [8930.693588789, 8902.132816716, 8877.254465461, 8855.095904154, 8835.402527951, 8816.99076143, 8799.611013471, 8782.171242747, 8767.422177564, 8750.683101494]
CHFBiasi = [4868.190469249, 4863.216092019, 4857.466714621, 4851.372467595, 4845.046740207, 4838.534306776, 4831.846994929, 4825.145122596, 4818.329380123, 4813.463682783]
CHFOKB = [12237.172162431, 12179.833238532, 12124.548260665, 12071.946824535, 12021.248231219, 11971.760170627, 11923.12579644, 11874.935545858, 11829.25791467, 11788.002744576]
CHFBowring = [8270.30333, 8250.821268, 8233.677158, 8218.56515, 8204.861225, 8192.029454, 8179.860549, 8167.001861, 8156.619242, 8145.728677]

CHFA = [CHFW3[i]*K3_B[i] for i in range(10)]
CHFB = [CHFBiasi[i]*K3_B[i] for i in range(10)]
CHFC = [CHFOKB[i]*K3_B[i] for i in range(10)]
CHFD = [CHFBowring[i]*K3_B[i] for i in range(10)]

Sum7 = 0
Sum8 = 0
Sum9 = 0
Sum10 = 0
for i in range(10):
    Sum7 = Sum7 + n.abs(CHF_BI[i] - CHFA[i])
    Sum8 = Sum8 + n.abs(CHF_BI[i] - CHFB[i])
    Sum9 = Sum9 + n.abs(CHF_BI[i] - CHFC[i])
    Sum10 = Sum10 + n.abs(CHF_BI[i] - CHFD[i])
ER7 = Sum7*100/10/av
ER8 = Sum8*100/10/av
ER9 = Sum9*100/10/av
ER10 = Sum10*100/10/av

# Normalized Root Mean Square Error (NRMSE)

Sum11 = 0
Sum12 = 0
Sum13 = 0
Sum14 = 0
for i in range(10):
    Sum11 = Sum11 + n.abs(CHF_BI[i] - CHFA[i])**2
    Sum12 = Sum12 + n.abs(CHF_BI[i] - CHFB[i])**2
    Sum13 = Sum13 + n.abs(CHF_BI[i] - CHFC[i])**2
    Sum14 = Sum14 + n.abs(CHF_BI[i] - CHFD[i])**2
ER11 = n.sqrt(Sum11/10)*100/av
ER12 = n.sqrt(Sum12/10)*100/av
ER13 = n.sqrt(Sum13/10)*100/av
ER14 = n.sqrt(Sum14/10)*100/av

CHF1 = n.array(CHF_GI2)
CHF2 = n.array(CHF_BI)
CHF3 = n.array(CHFA)
CHF4 = n.array(CHFB)
CHF5 = n.array(CHFC)
CHF6 = n.array(CHFD)

p.plot(N, CHF1, "--sb", label = "Groeneveld LUT. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER2, ER5))
p.plot(N, CHF2, "-.Dg", label = "Bobkov LUT")
p.plot(N, CHF3, "--*r", label = "W-3 Correlation. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER7, ER11))
p.plot(N, CHF4, "-.pc", label = "Biasi Correlation. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER8, ER12))
p.plot(N, CHF5, "--om", label = "OKB-Gidropress Correlation. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER9, ER13))
p.plot(N, CHF6, "-.^k", label = "Bowring Correlation. NMAE = %.2f%%. NRMSE = %.2f%%" % (ER10, ER14))

p.grid("on")
p.xlabel("Normalized length ($z/H$)")
p.ylabel("CHF ($kW/m^{2}$)")
p.legend()
p.show()
