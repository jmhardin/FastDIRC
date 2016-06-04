import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

theta_spacing=0.5-.000001
phi_spacing=1.0-.000001

theta_min=0
theta_max=12
phi_min=0
phi_max=90

filename="res_fullscan_fine_noceff.csv"
theta, phi, pion_y, kaon_y, kde, lin, lin_cal, lut = np.loadtxt(filename, delimiter=' ', usecols=(2, 3, 5, 7, 8, 9, 10, 11), unpack=True)

nphi = np.unique(phi).size
ntheta = np.unique(theta).size
KDE = np.zeros((ntheta,nphi)) 
LIN = np.zeros((ntheta,nphi)) 
LIN_CAL = np.zeros((ntheta,nphi)) 
LUT = np.zeros((ntheta,nphi)) 
LIN_CAL_IMP = np.zeros((ntheta,nphi)) 
LIN_CAL_V_KDE = np.zeros((ntheta,nphi)) 
LUT_V_KDE = np.zeros((ntheta,nphi)) 
LUT_V_LIN_CAL = np.zeros((ntheta,nphi)) 
FUDGE_PION = np.zeros((ntheta,nphi)) 
FUDGE_KAON = np.zeros((ntheta,nphi)) 

assert phi.shape == theta.shape == kde.shape
for i in range(len(phi)):
    KDE[theta[i]/theta_spacing][phi[i]/phi_spacing] = kde[i] 
    LIN[theta[i]/theta_spacing][phi[i]/phi_spacing] = lin[i] 
    LIN_CAL[theta[i]/theta_spacing][phi[i]/phi_spacing] = lin_cal[i] 
    LUT[theta[i]/theta_spacing][phi[i]/phi_spacing] = lut[i] 
    LIN_CAL_IMP[theta[i]/theta_spacing][phi[i]/phi_spacing] = (lin_cal[i]-lin[i])/lin[i]
    LIN_CAL_V_KDE[theta[i]/theta_spacing][phi[i]/phi_spacing] = (lin_cal[i]-kde[i])/kde[i]
    LUT_V_KDE[theta[i]/theta_spacing][phi[i]/phi_spacing] = (lut[i]-kde[i])/kde[i]
    LUT_V_LIN_CAL[theta[i]/theta_spacing][phi[i]/phi_spacing] = (lin_cal[i]-lut[i])/lin_cal[i]
    FUDGE_PION[theta[i]/theta_spacing][phi[i]/phi_spacing] = pion_y[i]
    FUDGE_KAON[theta[i]/theta_spacing][phi[i]/phi_spacing] = kaon_y[i]

phi_ticks = np.unique(phi)
theta_ticks = np.unique(theta)

clim_low = 1
clim_high = 3

phi_label="Phi (deg)"
theta_label="Theta (deg)"
res_label="Resolution (mrad)"
improvement_label="Percent Improvement"
worse_label="Percent Worse"
fudge_label="Y correction (mm)"


plt.pcolor(phi_ticks,theta_ticks,KDE,cmap=plt.cm.afmhot)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(res_label)
plt.savefig('kde_map.pdf')
plt.clf()
plt.pcolor(phi_ticks,theta_ticks,LIN,cmap=plt.cm.afmhot)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(res_label)
plt.savefig('lin_map.pdf')
plt.clf()
plt.pcolor(phi_ticks,theta_ticks,LIN_CAL,cmap=plt.cm.afmhot)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(res_label)
plt.savefig('lin_cal_map.pdf')
plt.clf()
plt.pcolor(phi_ticks,theta_ticks,LUT,cmap=plt.cm.afmhot)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(res_label)
plt.savefig('lut_map.pdf')
plt.clf()

clim_low = -1
clim_high = 0

plt.pcolor(phi_ticks,theta_ticks,LIN_CAL_IMP,cmap=plt.cm.afmhot)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(improvement_label)
plt.savefig('lin_cal_imp_map.pdf')
plt.clf()

clim_low = 0
clim_high = 1

plt.pcolor(phi_ticks,theta_ticks,LIN_CAL_V_KDE,cmap=plt.cm.afmhot)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(worse_label)
plt.savefig('lin_cal_v_kde_map.pdf')
plt.clf()

plt.pcolor(phi_ticks,theta_ticks,LUT_V_KDE,cmap=plt.cm.afmhot)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(worse_label)
plt.savefig('lut_v_kde_map.pdf')
plt.clf()

clim_low = -1
clim_high = 1

plt.pcolor(phi_ticks,theta_ticks,LUT_V_LIN_CAL,cmap=plt.cm.seismic)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(worse_label)
plt.savefig('lut_v_lin_cal_map.pdf')
plt.clf()

clim_low = -15
clim_high = 15

plt.pcolor(phi_ticks,theta_ticks,FUDGE_PION,cmap=plt.cm.seismic)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(fudge_label)
plt.savefig('pion_fudge_map.pdf')
plt.clf()

clim_low = -15
clim_high = 15

plt.pcolor(phi_ticks,theta_ticks,FUDGE_KAON,cmap=plt.cm.seismic)
plt.clim(clim_low,clim_high)
plt.axis([phi_min,phi_max,theta_min,theta_max])
plt.xlabel(phi_label)
plt.ylabel(theta_label)
cbar = plt.colorbar()
cbar.ax.set_ylabel(fudge_label)
plt.savefig('kaon_fudge_map.pdf')
plt.clf()

