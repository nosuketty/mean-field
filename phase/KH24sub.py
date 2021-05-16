import numpy as np
import matplotlib.pyplot as plt

texparams = {'ps.useafm': True, 'pdf.use14corefonts': True, 'pdf.fonttype': 42,
             'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsmath} \usepackage{txfonts} \usepackage{bm}'}
plt.rcParams.update(texparams)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

w = np.loadtxt('24sub.txt')
g = np.loadtxt('afstar24sub.txt')
f = np.loadtxt('fmstar24sub.txt')
v = np.loadtxt('vortex24sub_.txt')
n = np.loadtxt('fmneel.txt')
av = np.loadtxt('afvortex24sub_.txt')

plt.figure(figsize=(15, 8), dpi=50)


plt.plot(w[:,0], w[:,1], lw=0.5, ls='dashed', color='black')
plt.plot(g[:,0], g[:,1], lw=0.5, color='black')
plt.plot(f[:,0], f[:,1], lw=0.5, color='black')
plt.plot(v[:,0], v[:,1], lw=0.5, color='black')
plt.plot(n[:,0], n[:,1], lw=0.5, color='black')
plt.plot(av[:,0], av[:,1], lw=0.5, color='black')


plt.fill_between([-np.pi-0.03,np.pi+0.03],[9,9],facecolor='white',alpha=10)  #ffm
plt.fill_between(w[:,0],w[:,1],facecolor='#FFCCCC',alpha=10)  #neel
plt.fill_between(w[0:41,0],w[0:41,1], facecolor='#F2C300',alpha=10)  #vortex
plt.fill_between(v[:,0],v[:,1], facecolor='#8EE0FB',alpha=10)  #fmstar
plt.fill_between(w[39:50,0],w[39:50,1], facecolor='#8EE0FB',alpha=10)  #fmstar
plt.fill_between(n[:,0],n[:,1], facecolor='#8EE0FB',alpha=10)  #fmstar
plt.fill_between(v[0:57,0],v[0:57,1], facecolor='#C7F8DA',alpha=10)  #canted stripy
#plt.fill_between(v[29:,0],v[29:,1], facecolor='#FFCCCC',alpha=10)  #neel
plt.fill_between(f[:,0],f[:,1], facecolor='#C7F8DA',alpha=10)  #canted stripy
plt.fill_between(w[157:193,0],w[157:193,1], facecolor='#D6FFCC',alpha=10)  #afvortex
plt.fill_between(w[191:,0],w[191:,1], facecolor='#D6CCFF',alpha=10)  #canted zigzag
plt.fill_between(av[:,0],av[:,1], facecolor='#D6CCFF',alpha=10)  #canted zigzag
plt.fill_between(g[:,0],g[:,1], facecolor='#F5FFCC',alpha=10) #afstar

st = np.array([[-1.5708,0.03],[-0.47085,0.03]])
ne = np.array([[-0.47085,0.03],[1.5708,0.03]])
zi = np.array([[1.57,0.03],[2.67035,0.03]])
fm1 = np.array([[-np.pi,0.03],[-np.pi/2,0.03]])
fm2 = np.array([[np.pi/2,0.03],[np.pi,0.03]])

plt.fill_between(fm1[:,0],fm1[:,1], facecolor='black',alpha=10, zorder=10) #fm
plt.fill_between(fm2[:,0],fm2[:,1], facecolor='black',alpha=10, zorder=10) #fm
plt.fill_between(st[:,0],st[:,1], facecolor='#31F07D',alpha=10, zorder=10) #stripy
plt.fill_between(ne[:,0],ne[:,1], facecolor='#F03030',alpha=10, zorder=10) #neel
plt.fill_between(zi[:,0],zi[:,1], facecolor='#5730F1',alpha=10, zorder=10) #zigzag








#plt.legend()
plt.xlabel(r'$\alpha$',fontsize=26)
plt.ylabel(r'$h$',fontsize=26)
plt.yticks([0,1,2,3,4],[0,1,2,3,4],fontsize=20)
plt.xticks([-np.pi,-np.pi/2, 0, np.pi/2,np.pi],[r'$-\pi$',r'$-\pi/2$',0,r'$\pi/2$',r'$\pi$'],fontsize=20)

plt.vlines(x=np.pi/2, ymin=0, ymax=2.0, color='black', lw=0.5)

plt.xlim([-3.14,3.14])
plt.ylim([0,4.5])
# plt.xlim([-np.pi/2-0.03, -np.pi/2+0.03])
# plt.ylim([0,0.05])
# plt.yticks(np.arange(0, 0.051, 0.01), np.arange(0, 0.051, 0.01))

#plt.fill_between(w[:,0]/np.pi,0,fs[:,1]*np.sqrt(3),facecolor='#C5E0B4',alpha=10)

plt.savefig('24sub.pdf',transparent=True, bbox_inches='tight')