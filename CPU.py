import numpy as n
import matplotlib.pyplot as p
import time as t

# V_PXG

V000 = 8718 # P = 14, x = -0.5, G = 4500
V001 = 8978 # P = 14, x = -0.5, G = 5000
V010 = 7543 # P = 14, x = -0.4, G = 4500
V011 = 7649 # P = 14, x = -0.4, G = 5000
V100 = 8015 # P = 16, x = -0.5, G = 4500
V101 = 8249 # P = 16, x = -0.5, G = 5000
V110 = 6902 # P = 16, x = -0.4, G = 4500
V111 = 7179 # P = 16, x = -0.4, G = 5000

G1 = 4500
G2 = 5000

x1 = -0.5
x2 = -0.4

P1 = 14
P2 = 16

x = [-0.428252, -0.427085, -0.425905, -0.424745, -0.423599, -0.422459, -0.42132, -0.420184, -0.419079, -0.418145]
G = [4513.528727783, 4504.715675237, 4498.613443827, 4494.3914841, 4491.473025599, 4489.419722781, 4488.081615326, 4486.881932781, 4486.974216054, 4483.859655599]
P = [15.5737, 15.5653, 15.5569, 15.5485, 15.54, 15.5316, 15.5232, 15.5148, 15.5063, 15.4979]


G_x_p = []
G_p_x = []
p_G_x = []
p_x_G = []
x_G_p = []
x_p_G = []
Bourke = []

s = 10**6 # Number of computational cycles

for j in range(s):
#-------------------------------------------------------------------------------
# G then x then P
#-------------------------------------------------------------------------------

    t0 = t.time()

# Interpolations at P1

# Interpolation at x1


    A01 = []
    for i in range(10):
        a = ((G[i] - G1) * (V001 - V000) / (G2 - G1)) + V000
        A01.append(a)

# Interpolation at x2
    A02 = []
    for i in range(10):
        a = ((G[i] - G1) * (V011 - V010) / (G2 - G1)) + V010
        A02.append(a)

# Interpolation at G
    A03 = []
    for i in range(10):
        a = ((x[i] - x1) * (A02[i] - A01[i]) / (x2 - x1)) + A01[i]
        A03.append(a)
#-------------------------------------------------------------------------------
# Interpolations at P2

# Interpolation at x1
    A11 = []
    for i in range(10):
        a = ((G[i] - G1) * (V101 - V100) / (G2 - G1)) + V100
        A11.append(a)

# Interpolation at x2
    A12 = []
    for i in range(10):
        a = ((G[i] - G1) * (V111 - V110) / (G2 - G1)) + V110
        A12.append(a)

# Interpolation at G
    A13 = []
    for i in range(10):
        a = ((x[i] - x1) * (A12[i] - A11[i]) / (x2 - x1)) + A11[i]
        A13.append(a)

# Final values of the CHF
    A = []
    for i in range(10):
        a = ((P[i] - P1) * (A13[i] - A03[i]) / (P2 - P1)) + A03[i]
        A.append(a)
#-------------------------------------------------------------------------------
# G then P then x
#-------------------------------------------------------------------------------
    t1 = t.time()

    dt1 = t1 - t0

# Interpolations at x1

# Interpolation at P1
    B01 = []
    for i in range(10):
        a = ((G[i] - G1) * (V001 - V000) / (G2 - G1)) + V000
        B01.append(a)

# Interpolation at P2
    B02 = []
    for i in range(10):
        a = ((G[i] - G1) * (V101 - V100) / (G2 - G1)) + V100
        B02.append(a)

# Interpolation at G w.r.t. P
    B03 = []
    for i in range(10):
        a = ((P[i] - P1) * (B02[i] - B01[i]) / (P2 - P1)) + B01[i]
        B03.append(a)
#-------------------------------------------------------------------------------
# Interpolations at x2

# Interpolation at P1
    B11 = []
    for i in range(10):
        a = ((G[i] - G1) * (V011 - V010) / (G2 - G1)) + V010
        B11.append(a)

# Interpolation at P2
    B12 = []
    for i in range(10):
        a = ((G[i] - G1) * (V111 - V110) / (G2 - G1)) + V110
        B12.append(a)

# Interpolation at G
    B13 = []
    for i in range(10):
        a = ((P[i] - P1) * (B12[i] - B11[i]) / (P2 - P1)) + B11[i]
        B13.append(a)

# Final values of the CHF
    B = []
    for i in range(10):
        a = ((x[i] - x1) * (B13[i] - B03[i]) / (x2 - x1)) + B03[i]
        B.append(a)
#-------------------------------------------------------------------------------
# P then G then x
#-------------------------------------------------------------------------------
    t2 = t.time()

    dt2 = t2 - t1

# Interpolations at x1

# Interpolation at G1
    C01 = []
    for i in range(10):
        a = ((P[i] - P1) * (V100 - V000) / (P2 - P1)) + V000
        C01.append(a)

# Interpolation at G2
    C02 = []
    for i in range(10):
        a = ((P[i] - P1) * (V101 - V001) / (P2 - P1)) + V001
        C02.append(a)

# Interpolation at P
    C03 = []
    for i in range(10):
        a = ((G[i] - G1) * (C02[i] - C01[i]) / (G2 - G1)) + C01[i]
        C03.append(a)
#-------------------------------------------------------------------------------
# Interpolations at x2

# Interpolation at G1
    C11 = []
    for i in range(10):
        a = ((P[i] - P1) * (V110 - V010) / (P2 - P1)) + V010
        C11.append(a)

# Interpolation at G2
    C12 = []
    for i in range(10):
        a = ((P[i] - P1) * (V111 - V011) / (P2 - P1)) + V011
        C12.append(a)

# Interpolation at P
    C13 = []
    for i in range(10):
        a = ((G[i] - G1) * (C12[i] - C11[i]) / (G2 - G1)) + C11[i]
        C13.append(a)

# Final values of the CHF
    C = []
    for i in range(10):
        a = ((x[i] - x1) * (C13[i] - C03[i]) / (x2 - x1)) + C03[i]
        C.append(a)
#-------------------------------------------------------------------------------
# P then x then G
#-------------------------------------------------------------------------------
    t3 = t.time()

    dt3 = t3 - t2

# Interpolations at G1

# Interpolation at x1
    D01 = []
    for i in range(10):
        a = ((P[i] - P1) * (V100 - V000) / (P2 - P1)) + V000
        D01.append(a)

# Interpolation at x2
    D02 = []
    for i in range(10):
        a = ((P[i] - P1) * (V110 - V010) / (P2 - P1)) + V010
        D02.append(a)

# Interpolation at P
    D03 = []
    for i in range(10):
        a = ((x[i] - x1) * (D02[i] - D01[i]) / (x2 - x1)) + D01[i]
        D03.append(a)
#-------------------------------------------------------------------------------
# Interpolations at G2

# Interpolation at x1
    D11 = []
    for i in range(10):
        a = ((P[i] - P1) * (V101 - V001) / (P2 - P1)) + V001
        D11.append(a)

# Interpolation at x2
    D12 = []
    for i in range(10):
        a = ((P[i] - P1) * (V111 - V011) / (P2 - P1)) + V011
        D12.append(a)

# Interpolation at P
    D13 = []
    for i in range(10):
        a = ((x[i] - x1) * (D12[i] - D11[i]) / (x2 - x1)) + D11[i]
        D13.append(a)

# Final values of the CHF
    D = []
    for i in range(10):
        a = ((G[i] - G1) * (D13[i] - D03[i]) / (G2 - G1)) + D03[i]
        D.append(a)
#-------------------------------------------------------------------------------
# x then G then P
#-------------------------------------------------------------------------------
    t4 = t.time()

    dt4 = t4 - t3

# Interpolations at P1

# Interpolation at G1
    E01 = []
    for i in range(10):
        a = ((x[i] - x1) * (V010 - V000) / (x2 - x1)) + V000
        E01.append(a)

# Interpolation at G2
    E02 = []
    for i in range(10):
        a = ((x[i] - x1) * (V011 - V001) / (x2 - x1)) + V001
        E02.append(a)

# Interpolation at x
    E03 = []
    for i in range(10):
        a = ((G[i] - G1) * (E02[i] - E01[i]) / (G2 - G1)) + E01[i]
        E03.append(a)
#-------------------------------------------------------------------------------
# Interpolations at P2

# Interpolation at G1
    E11 = []
    for i in range(10):
        a = ((x[i] - x1) * (V110 - V100) / (x2 - x1)) + V100
        E11.append(a)

# Interpolation at G2
    E12 = []
    for i in range(10):
        a = ((x[i] - x1) * (V111 - V101) / (x2 - x1)) + V101
        E12.append(a)

# Interpolation at x
    E13 = []
    for i in range(10):
        a = ((G[i] - G1) * (E12[i] - E11[i]) / (G2 - G1)) + E11[i]
        E13.append(a)

# Final values of the CHF
    E = []
    for i in range(10):
        a = ((P[i] - P1) * (E13[i] - E03[i]) / (P2 - P1)) + E03[i]
        E.append(a)
#-------------------------------------------------------------------------------
# x then P then G
#-------------------------------------------------------------------------------
    t5 = t.time()

    dt5 = t5 - t4

# Interpolations at G1

# Interpolation at P1
    F01 = []
    for i in range(10):
        a = ((x[i] - x1) * (V010 - V000) / (x2 - x1)) + V000
        F01.append(a)

# Interpolation at P2
    F02 = []
    for i in range(10):
        a = ((x[i] - x1) * (V110 - V100) / (x2 - x1)) + V100
        F02.append(a)

# Interpolation at x
    F03 = []
    for i in range(10):
        a = ((P[i] - P1) * (F02[i] - F01[i]) / (P2 - P1)) + F01[i]
        F03.append(a)
#-------------------------------------------------------------------------------
# Interpolations at G2

# Interpolation at P1
    F11 = []
    for i in range(10):
        a = ((x[i] - x1) * (V011 - V001) / (x2 - x1)) + V001
        F11.append(a)

# Interpolation at P2
    F12 = []
    for i in range(10):
        a = ((x[i] - x1) * (V111 - V101) / (x2 - x1)) + V101
        F12.append(a)

# Interpolation at x
    F13 = []
    for i in range(10):
        a = ((P[i] - P1) * (F12[i] - F11[i]) / (P2 - P1)) + F11[i]
        F13.append(a)

# Final values of the CHF
    F = []
    for i in range(10):
        a = ((G[i] - G1) * (F13[i] - F03[i]) / (G2 - G1)) + F03[i]
        F.append(a)
#-------------------------------------------------------------------------------
# A Trilinear Interpolation Algorithm by Paul Bourke
# Link: http://paulbourke.net/miscellaneous/interpolation/

# a = scaled pressure
# b = scaled equilibrium quality
# c = scaled mass flux
# V_xyz = V_abc = V_PXG

# Scaling pressure

    t6 = t.time()

    dt6 = t6 - t5

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

    CHF = []
    for i in range(10):
        z = V000*(1 - a[i])*(1 - b[i])*(1 - c[i]) + V100*a[i]*(1 - b[i])*(1 - c[i]) + V010*(1 - a[i])*b[i]*(1 - c[i]) + V001*(1 - a[i])*(1 - b[i])*c[i] + V101*a[i]*(1 - b[i])*c[i] + V011*(1 - a[i])*b[i]*c[i] + V110*a[i]*b[i]*(1 - c[i]) + V111*a[i]*b[i]*c[i]
        CHF.append(z)

    t7 = t.time()

    dt7 = t7 - t6

    G_x_p.append(dt1*10**6)
    G_p_x.append(dt2*10**6)
    p_G_x.append(dt3*10**6)
    p_x_G.append(dt4*10**6)
    x_G_p.append(dt5*10**6)
    x_p_G.append(dt6*10**6)
    Bourke.append(dt7*10**6)

#===============================================================================

A = n.array(G_x_p)
B = n.array(G_p_x)
C = n.array(p_G_x)
D = n.array(p_x_G)
E = n.array(x_G_p)
F = n.array(x_p_G)
G = n.array(Bourke)

means = [n.mean(A), n.mean(B), n.mean(C), n.mean(D), n.mean(E), n.mean(F), n.mean(G)]
standard_deviations = [n.std(A), n.std(B), n.std(C), n.std(D), n.std(E), n.std(F), n.std(G)]

print("Means:", means)
print("Standard Deviations:", standard_deviations)
"""

N1 = n.array([i+1 for i in range(s)])

p.plot(N1, A, "--m", label = "$G-x_{e}-p$")
p.plot(N1, B, ":m", label = "$G-p-x_{e}$")
p.plot(N1, C, "--g", label = "$p-G-x_{e}$")
p.plot(N1, D, ":g", label = "$p-x_{e}-G$")
p.plot(N1, E, "--c", label = "$x_{e}-G-p$")
p.plot(N1, F, ":c", label = "$x_{e}-p-G$")
p.plot(N1, G, "--r", label = "Bourke Algorithm")

p.grid("on")
p.xlabel("Number of computational cycles")
p.ylabel("CPU time ($\mu$s)")
p.legend()
p.show()

# Inclusion of K1 and K4

K1 = (8/10.5)**(0.425)

K1CHF = []
for i in range(10):
    a = K1*CHF[i]
    K1CHF.append(a)

alpha = [1.792456951, 1.794026504, 1.796968338, 1.799242769, 1.801253564, 1.803205343, 1.805172303, 1.807145907, 1.809112296, 1.811060236]

K4 = []
for i in range(10):
    L = 0.353*(i+1)
    j = n.exp((0.0105/L)*n.exp(2*alpha[i]))
    K4.append(j)

A1 = []
B1 = []
C1 = []
D1 = []
E1 = []
F1 = []
CHFPB = []
for i in range(10):
    a = K1*K4[i]*A[i]
    b = K1*K4[i]*B[i]
    c = K1*K4[i]*C[i]
    d = K1*K4[i]*D[i]
    e = K1*K4[i]*E[i]
    f = K1*K4[i]*F[i]
    g = K1*K4[i]*CHF[i]
    A1.append(a)
    B1.append(b)
    C1.append(c)
    D1.append(d)
    E1.append(e)
    F1.append(f)
    CHFPB.append(g)

CHFW3 = [8930.693588789, 8900.219989696, 8873.663343745, 8849.935960247, 8828.538109008, 8808.441650409, 8789.381493, 8770.582505548, 8754.124545518, 8741.098010377]
CHFBiasi = [4868.190469249, 4863.216092019, 4857.466714621, 4851.372467595, 4845.046740207, 4838.534306776, 4831.846994929, 4825.145122596, 4818.329380123, 4812.672639425]
CHFVVER = [12237.172162431, 12179.833238532, 12124.548260665, 12071.946824535, 12021.248231219, 11971.760170627, 11923.12579644, 11874.935545858, 11829.25791467, 11788.002744576]
CHFW3ZIS = [7632.124148745, 7608.522462743, 7588.063418671, 7570.19255973, 7554.157479847, 7539.214876132, 7525.151653591, 7511.299619408, 7499.48019026, 7486.220542109]

N = n.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
AA = n.array(A1)
BB = n.array(B1)
CC = n.array(C1)
DD = n.array(D1)
EE = n.array(E1)
FF = n.array(F1)
CHF1 = n.array(CHFPB)
CHF2 = n.array(CHFW3)
CHF3 = n.array(CHFBiasi)
CHF4 = n.array(CHFVVER)
H1 = n.array(K1CHF)
H2 = n.array(CHF)
H3 = n.array(CHFW3ZIS)

#===============================================================================

# RESULTS

# print("Standard Algorithm:", dt1*10**6, dt2*10**6, dt3*10**6, dt4*10**6, dt5*10**6, dt6*10**6)
# print("Bourke Algorithm:", dt7*10**6)

# p.plot(N, AA, "--r", label = "$G-x_{e}-p$ Curve")
# p.plot(N, BB, "--r", label = "$G-p-x_{e}$ Curve")
# p.plot(N, CC, "--g", label = "$p-G-x_{e}$ Curve")
# p.plot(N, DD, "--g", label = "$p-x_{e}-G$ Curve")
# p.plot(N, EE, "--c", label = "$x_{e}-G-p$ Curve")
# p.plot(N, FF, "--c", label = "$x_{e}-p-G$ Curve")

p.plot(N, CHF1, "-k", label = "$K_1 \ K_4 \ {CHF}_{table}$")
p.plot(N, H1, "--g", label = "$K_1 \ {CHF}_{table}$")
p.plot(N, H2, "--r", label = "${CHF}_{table}$")
p.plot(N, [3000 for i in range(10)], " ")

p.grid("on")
p.xlabel("Normalized length ($z/H$)")
p.ylabel("CHF ($kW/m^{2}$)")
p.legend()
p.show()

#===============================================================================

p.plot(N, CHF1, "-k", label = "$K_1 \ K_4 \ {CHF}_{table}$")
p.plot(N, CHF2, "--g", label = "W-3 Correlation")
p.plot(N, H3, "--m", label = "W-3 Correlation with $\Delta h_{sub}$ = 0")
p.plot(N, CHF3, "--r", label = "Biasi et al. Correlation")
p.plot(N, CHF4, "--c", label = "OKB-Gidorpress Correlation")

p.plot(N, [3000 for i in range(10)], " ")

p.grid("on")
p.xlabel("Normalized length ($z/H$)")
p.ylabel("CHF ($kW/m^{2}$)")
p.legend()
p.show()
"""
