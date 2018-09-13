"""Finite Differences method
   for the Heat Equation (Laplace 2D):
-----------------------------------------------------------------------------
----------               u_t = alpha * (uxx + u_yy)                ----------
-----------------------------------------------------------------------------
Boundary conditions assumed:
                                    du/dx=0
                                --------------
                                |            |
                                |            |
                        du/dx=0 |   random   | du/dx=0
                                |   (t=0)    |
                                |            |
                                --------------
                                    du/dx=0
-----------------------------------------------------------------------------
Observations: code with 70 columns.
"""
__version__ = 1.0
__author__ = """Daniel Coelho (danielcoelho.uerj@gmail.com)"""
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
import numpy as np
# ----------------------------------------------------------------------------
L = H = .7
T_initial = 0  # deg C, initial temp of wall
T1s = 0  # deg C, surface 1 temp
T2s = 50  # deg C, surface 2 temp
alpha = 2.8e-4  # thermal dif
beta = 5e-3  # thermal dif
tauinv = 50.
k = -.005
t_final = 10.  # input("Final time? ")              #s, simulation time
dt = 5e-2  # input("Time Step? ")                    #s, fixed time step
n = int(t_final / dt)
print
n
# ----------------------------------------------------------------------------
# np.random.seed(0)
"""Testando matrizes"""
Ts = 0.
Ts2 = 0.
T0 = 0.
Q = 0.
nx = ny = 64
dx = L / (nx - 2)  # m, node thickness
dy = dx  # H/(ny-2)                        #m, node thickness
h = dx
print('dx: ' + str(dx))
print('dy: ' + str(dy))
print('dt: ' + str(dt))
x = np.linspace(-0.5 * dx, L + 0.5 * dx, nx, dtype='float')
y = np.linspace(-0.5 * dy, H + 0.5 * dy, ny, dtype='float')
t = np.arange(0.0, t_final, dt, dtype='float')
# un=Ts*np.ones((int(ny)))
# Tbulk=T0*np.ones((int(ny)))
un = np.random.rand(nx, ny)
# un=10*un
# ----------------------------------------------------------------------------
"""Testando matrizes"""
vs = 0.
vs2 = 0.
v0 = 0.
vn = vs * np.ones((int(ny)))
vbulk = v0 * np.ones((int(ny)))
vn = np.random.rand(nx, ny)
# for Z in (un,vn):
#     Z[0,:]=Z[1,:]
#     Z[-1,:]=Z[-2,:]
#     Z[:,0]=Z[:,1]
#     Z[:,-1]=Z[:,-2]
un_i = un.copy()
vn_i = vn.copy()
U0 = un.copy()
V0 = vn.copy()
print('len(un): ' + str(len(un)))
print('len(vn): ' + str(len(vn)))
print('Maior valor para a especie un: ' +
      str(un.max()))
print('Menor valor para a especie un: ' +
      str(un.min()))
print('Maior valor para a especie vn: ' +
      str(vn.max()))
print('Menor valor para a especie vn: ' +
      str(vn.min()))
# ----------------------------------------------------------------------------
A = 0.5 * alpha * dt / h ** 2
B = 2 * A
# ----------------------------------------------------------------------------
A1 = 0.5 * tauinv * beta * dt / h ** 2
B1 = 2 * A1
# ----------------------------------------------------------------------------
mfig = [0, int(1. / 7 * n), int(2. / 7 * n), int(3. / 7 * n),
        int(4. / 7 * n), int(5. / 7 * n), int(6. / 7 * n),
        int(7. / 7 * n) - 1]
fignum = 0
unovo1 = np.zeros((ny, nx), dtype='float')
unovo2 = np.zeros((ny, nx), dtype='float')
vnovo1 = np.zeros((ny, nx), dtype='float')
vnovo2 = np.zeros((ny, nx), dtype='float')
bstep1 = np.zeros((nx), dtype='float')
bstep2 = np.zeros((ny), dtype='float')
o = 0
normL1 = []
e = 10 ** (-5)
u = 0
for o in range(0, len(t)):
    # ----------------------------------------------------------------------------
    #                                 Species u                                 -
    # ----------------------------------------------------------------------------
    # ----------              First Splitting Equation                  ----------
    K = np.zeros((ny, nx), dtype='float')
    LambdaY = np.zeros((ny - 1, nx - 1), dtype='float')
    LambdaYun = np.zeros((ny - 1, nx - 1), dtype='float')
    bstep1 = np.zeros((nx), dtype='float')
    K[0, 1] = -A
    K[0, nx - 1] = -A
    K[nx - 1, 0] = -A
    K[nx - 1, nx - 2] = -A
    for i in range(1, ny - 1):
        for j in range(1, nx - 1):
            K[0, 0] = 1 + B + 0.5 * dt * (un[i, j] ** 2)
            K[nx - 1, nx - 1] = K[0, 0]
            K[j, j - 1] = -A
            K[j, j] = 1 + B + 0.5 * dt * (un[i, j] ** 2)
            K[j, j + 1] = -A
            fn = un[i, j] - vn[i, j] + k
            LambdaY[i, j] = -B - 0.5 * dt * (un[i, j] ** 2)
            LambdaYun[i, j] = (
                    A * un[i - 1, j] +
                    LambdaY[i, j] * un[i, j] +
                    A * un[i + 1, j]
            )
            bstep1[j] = un[i, j] + LambdaYun[i, j] + dt * fn
        #         bstep1[0]=un[i,0]+LambdaYun[i,0]+dt*fn
        #         bstep1[nx-1]=un[i,nx-1]+LambdaYun[i,nx-1]+dt*fn
        unovo1[i, :] = np.linalg.solve(K, bstep1)
    #         unovo2=unovo1.copy()
    # ----------              Second Splitting Equation                  ----------
    bstep2 = np.zeros((ny), dtype='float')
    for j in range(1, nx - 1):
        for i in range(1, ny - 1):
            K[i, i - 1] = -A
            K[i, i] = 1 - LambdaY[i, j]
            K[i, i + 1] = -A
            bstep2[i] = unovo1[i, j] - LambdaYun[i, j]
        #         bstep2[0]=unovo1[0,j]-LambdaYun[0,j]
        #         bstep2[nx-1]=unovo1[nx-1,j]-LambdaYun[nx-1,j]
        unovo2[:, j] = np.linalg.solve(K, bstep2)
    unovo2[0, 0] = 0.5 * (unovo2[1, 0] + unovo2[0, 1])
    unovo2[0, nx - 1] = 0.5 * (unovo2[1, nx - 1] + unovo2[0, nx - 2])
    unovo2[ny - 1, 0] = 0.5 * (unovo2[ny - 2, 0] + unovo2[ny - 1, 1])
    unovo2[ny - 1, ny - 1] = 0.5 * (unovo2[ny - 1, ny - 2] + unovo2[ny - 2, ny - 1])
    # ----------------------------------------------------------------------------
    #                                 Species v                                 -
    # ----------------------------------------------------------------------------
    # ----------              First Splitting Equation                  ----------
    K = np.zeros((ny, nx), dtype='float')
    LambdaY = np.zeros((ny, nx), dtype='float')
    LambdaYvn = np.zeros((ny, nx), dtype='float')
    bstep1 = np.zeros((nx), dtype='float')
    K[0, 0] = 1 + B1 + 0.5 * dt * tauinv
    K[0, 1] = -A1
    K[0, nx - 1] = -A1
    K[nx - 1, 0] = -A1
    K[nx - 1, nx - 1] = K[0, 0]
    K[nx - 1, nx - 2] = -A1
    for i in range(1, ny - 1):
        for j in range(1, nx - 1):
            K[j, j - 1] = -A1
            K[j, j] = 1 + B1 + 0.5 * dt * tauinv
            K[j, j + 1] = -A1
            fn = tauinv * un[i, j]
            LambdaY[i, j] = -B1 - 0.5 * dt * tauinv
            LambdaYvn[i, j] = A1 * vn[i - 1, j] + LambdaY[i, j] * vn[i, j] + A1 * vn[i + 1, j]
            bstep1[j] = vn[i, j] + LambdaYvn[i, j] + dt * fn
        #         bstep1[0]=vn[i,0]+LambdaYvn[i,0]+dt*fn
        #         bstep1[nx-1]=vn[i,nx-1]+LambdaYvn[i,nx-1]+dt*fn
        vnovo1[i, ::] = np.linalg.solve(K, bstep1)
    # ----------              Second Splitting Equation                  ----------
    bstep2 = np.zeros((ny), dtype='float')
    for j in range(1, nx - 1):
        for i in range(1, ny - 1):
            K[i, i - 1] = -A1
            K[i, i] = 1 - LambdaY[i, j]
            K[i, i + 1] = -A1
            bstep2[i] = vnovo1[i, j] - LambdaYvn[i, j]
            # -A*vn[i-1,j]+(B)*vn[i,j]-A*vn[i+1,j]
        #         bstep2[0]=vnovo1[0,j]-LambdaYvn[0,j]
        #         bstep2[nx-1]=vnovo1[nx-1,j]-LambdaYvn[nx-1,j]
        vnovo2[::, j] = np.linalg.solve(K, bstep2)
    vnovo2[0, 0] = 0.5 * (vnovo2[1, 0] + vnovo2[0, 1])
    vnovo2[0, nx - 1] = 0.5 * (vnovo2[1, nx - 1] + vnovo2[0, nx - 2])
    vnovo2[ny - 1, 0] = 0.5 * (vnovo2[ny - 2, 0] + vnovo2[ny - 1, 1])
    vnovo2[ny - 1, ny - 1] = 0.5 * (vnovo2[ny - 1, ny - 2] + vnovo2[ny - 2, ny - 1])
    # ----------------------------------------------------------------------------
    #                               Boundaries Fix                              -
    # ----------------------------------------------------------------------------
    #     for Z in (unovo2,vnovo2):
    # #         Z[0,:]=Z[1,:]
    # #         Z[-1,:]=Z[-2,:]
    #         Z[:,0]=Z[:,1]
    #         Z[:,-1]=Z[:,-2]
    # ----------------------------------------------------------------------------
    #                                 Norm L1                                   -
    # ----------------------------------------------------------------------------
    L1 = (1 / dt) * (abs(unovo2 - un).sum() / (abs(unovo2).sum()))
    normL1 = np.append(normL1, L1)
    # ----------------------------------------------------------------------------
    #                                 Results                                   -
    # ----------------------------------------------------------------------------
    if o in mfig:
        fignum += 1
        print(o, fignum)
        # if o == 0:
            # np.savetxt(str(dt) + '_' + str(o) + '.txt', un_i)
        # if o >= 1:
        #     #             np.savetxt('%d%.3dplot%d.csv'%(nx,dt,fignum),unovo2, delimiter=",")
        #     # np.savetxt(str(dt) + '_' + str(o) + '.txt', unovo2)
        #     T_ult = unovo2.copy()
    if L1 < e:
        print('-------------BREAK!!!-------------')
        break
    #     np.clip(unovo2,0,1,unovo2)
    #     np.clip(vnovo2,0,1,vnovo2)
    un = unovo2.copy()
    vn = vnovo2.copy()

print('Maior valor para a especie u: ' +
      str(unovo2.max()))
print('Menor valor para a especie u: ' +
      str(unovo2.min()))
print('Maior valor para a especie v: ' +
      str(vnovo2.max()))
print('Menor valor para a especie v: ' +
      str(vnovo2.min()))

