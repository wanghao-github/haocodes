from cmath import nan
from turtle import color, width
from unittest import result
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3d

j_data = np.loadtxt(r"C:\Users\wangh\Documents\HaoWang\Hao_code\不同U值YMnSn的J1J2J3\YMnSn_FM_AFM_123.txt")
U = j_data[:,0]
FM = j_data[:,1]
AFM1 = j_data[:,2]
AFM2 = j_data[:,3]
AFM3 = j_data[:,4]

co_matrix = np.array([[-1,-1,-2,1],[1,1,-2,1],[-1,1,2,1],[1,-1,2,1]])

j1j2j3E = np.zeros((27,4))
alpha_beta = np.zeros((27,2))
j2_j1 = np.zeros((27))
j3_j1 = np.zeros((27))
for u in np.arange(27):
    b = j_data[u,1:5]
    # print(b)
    result_J = np.linalg.solve(co_matrix,b.T)
    # print(result_J)
    j1j2j3E[u,:] = result_J
    j1 = j1j2j3E[u,0]
    j2 = j1j2j3E[u,1]
    j3 = j1j2j3E[u,2]
#     if j2*j3/(j1*j1)-j3/j2-j2/(4*j3)>=-1 and j2*j3/(j1*j1)-j3/j2-j2/(4*j3)<= 1:
#         alpha = -1*np.sign(j1*j3)*np.arccos(j2*j3/(j1*j1)-j3/j2-j2/(4*j3))
#     else:
#         alpha = 0
#     # beta = np.arccos(j1*j2/(8*j3*j3)-j2/(2*j1)-j1/(2*j2)) -alpha
#     if j3*j1/(j2*j2)-j1/(4*j3)-j3/j1 >=-1 and j3*j1/(j2*j2)-j1/(4*j3)-j3/j1 <=1:
#         beta = np.arccos(j3*j1/(j2*j2)-j1/(4*j3)-j3/j1)
#     else:
#         beta = 0
#     # print(alpha,beta)
#     alpha_beta[u,0] = alpha
#     alpha_beta[u,1] = beta
    j2_j1[u] = j2/j1
    j3_j1[u] = j3/j1

# np.savetxt(r"C:\Users\wangh\Documents\HaoWang\Hao_code\不同U值YMnSn的J1J2J3\YMnSn_FM_AFM_123.txt",a)
plt.figure(figsize=(20,18))

x=np.linspace(-3,0,500)
y=np.linspace(0,2,500)
def f1(x,y):
    return x + 2*y*(x+1)
def f2(x,y):
    return x - 2*y*(x+1)
def f3(x,y):
    return x - 2*y*(x-1)
def f4(x,y):
    return x + 2*y*(x-1)
X,Y=np.meshgrid(x,y)

plt.xlim(-1.2,0.2)
plt.ylim(-0.1,0.45)
plt.tick_params(labelsize=36,size=10,width=2)

ax = plt.gca()
ax.spines['left'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('right')
ax.spines['bottom'].set_linewidth(5)
ax.spines['right'].set_linewidth(5)
ax.spines['bottom'].set_position(('data',0))
ax.spines['right'].set_position(('data',0))

plt.contour(X,Y,f1(X,Y),0,linewidths=5)
plt.contour(X,Y,f2(X,Y),0,linewidths=5)
plt.contour(X,Y,f3(X,Y),0,linewidths=5)
plt.contour(X,Y,f4(X,Y),0,linewidths=5)
for i in np.arange(27):
    plt.scatter(j2_j1[i],j3_j1[i],s=600,color='red')
    plt.annotate(j_data[i,0],(j2_j1[i],j3_j1[i]),xytext=(j2_j1[i]+0.02,j3_j1[i]-0.002),color='red',fontsize=22)
    # plt.plot(j2_j1[i],j3_j1[i],'.',color='red')
    # plt.plot(j_data[i,0],(180/np.pi)*alpha_beta[i,0],'.',color='red')
    # plt.plot(j_data[i,0],(180/np.pi)*alpha_beta[i,1],'.',color='blue')
    # plt.plot(j_data[i,0],j1j2j3E[i,0],'.',color='red')
    # plt.plot(j_data[i,0],j1j2j3E[i,1],'.',color='blue')
    # plt.plot(j_data[i,0],j1j2j3E[i,2],'.',color='black')
    # plt.plot(j_data[i,0],j1j2j3E[i,3],'.',color='red')
plt.savefig(r'C:\Users\wangh\Desktop\YMnSn_U.png',dpi=600)
plt.show()