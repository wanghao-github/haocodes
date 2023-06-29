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

co_matrix = np.array([[1,1,2,1],[-1,-1,2,1],[1,-1,-2,1],[-1,1,-2,1]])

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

    a=j2*j3/(j1*j1)-j3/j2-j2/(4*j3)
    b=j3*j1/(j2*j2)-j1/(4*j3)-j3/j1
    if j2*j3/(j1*j1)-j3/j2-j2/(4*j3)>=-1 and j2*j3/(j1*j1)-j3/j2-j2/(4*j3)<= 1:
        alpha = -1*np.sign(j1*j3)*np.arccos(j2*j3/(j1*j1)-j3/j2-j2/(4*j3))
    else:
        alpha = 0
    # beta = np.arccos(j1*j2/(8*j3*j3)-j2/(2*j1)-j1/(2*j2)) -alpha
    if j3*j1/(j2*j2)-j1/(4*j3)-j3/j1 >=-1 and j3*j1/(j2*j2)-j1/(4*j3)-j3/j1 <=1:
        beta = np.arccos(j3*j1/(j2*j2)-j1/(4*j3)-j3/j1)
    else:
        beta = 0
    # print(a,b)
    # print(alpha,beta)
    alpha_beta[u,0] = alpha
    alpha_beta[u,1] = beta
    j2_j1[u] = j2/j1
    j3_j1[u] = j3/j1
np.savetxt(r"C:\Users\wangh\Desktop\allmightysumit3",j1j2j3E)
np.savetxt(r"C:\Users\wangh\Desktop\allmightysumit4",j_data)
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


np.savetxt(r"C:\Users\wangh\Desktop\kucha",f1(X,Y))

ax = plt.gca()

ax.spines['bottom'].set_linewidth(5)
ax.spines['right'].set_linewidth(5)
ax.spines['top'].set_linewidth(5)
ax.spines['left'].set_linewidth(5)


plt.contour(X,Y,f1(X,Y),0,linewidths=5)
plt.contour(X,Y,f2(X,Y),0,linewidths=5)
plt.contour(X,Y,f3(X,Y),0,linewidths=5)
plt.contour(X,Y,f4(X,Y),0,linewidths=5)
plt.axhline(y=0, ls='--',lw=5, c='red')
# plt.xlim(0,4.2)
# plt.ylim(-50,200)
plt.tick_params(labelsize=48,size=12,width=5)
for i in np.arange(27):
    # if i <=5:
        # plt.scatter(j_data[i,0],0,color='red',s=600)
        # plt.scatter(j_data[i,0],180,color='blue',s=600)
        # plt.tick_params(labelsize=48,size=12,width=5)
    plt.scatter(j2_j1[i],j3_j1[i],s=600,color='red')
    plt.annotate(j_data[i,0],(j2_j1[i],j3_j1[i]),xytext=(j2_j1[i]+0.02,j3_j1[i]-0.002),color='red',fontsize=22)
    # plt.plot(j2_j1[i],j3_j1[i],'.',color='red')
    # else:
    #     plt.scatter(j_data[i,0],(180/np.pi)*alpha_beta[i,0],color='red',s=600)
    #     plt.scatter(j_data[i,0],(180/np.pi)*alpha_beta[i,1],color='blue',s=600)
    #     plt.tick_params(labelsize=48,size=12,width=5)
    # plt.scatter(j_data[i,0],500*j1j2j3E[i,0],color='red',s=600)
    # plt.scatter(j_data[i,0],500*j1j2j3E[i,1],color='blue',s=600)
    # plt.scatter(j_data[i,0],500*j1j2j3E[i,2],color='black',s=600)
    # plt.scatter(j_data[i,0],j1j2j3E[i,3],color='red')
# plt.plot(j_data[i,0].tolist(),np.zeros((1,27)).tolist(),'.',color='red')


plt.savefig(r'C:\Users\wangh\Desktop\YMnSn_alpha_beta.png',dpi=600)
plt.show()