from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

name='WPMP'
structure=['st0']
namedist=['p_397C']
mn=0
mx=10
lev=15

for nd,struc in (zip(namedist,structure)):
######## WP and MP ##################
    #str0
    imgx=[]
    imgx.append(np.loadtxt('{nd}_wt_Glu0_{struc}.dat'.format(nd=nd,struc=struc)))
    imgx.append(np.loadtxt('{nd}_mut_Glu0_{struc}.dat'.format(nd=nd,struc=struc)))
    
    
    kernel = Gaussian2DKernel(x_stddev=1)
    X, Y = sorted(set(imgx[0][:,0])), sorted(set(imgx[0][:,1]))
    Z=[]
    #Glu0 wt and Glu0 mut
    for n in np.arange(2):
        print(n)
        Z.append(imgx[n][:,2].reshape(len(X), len(Y)))
        
    zz=[]
    #for z I make an avarage of the two strucs
    for n in np.arange(2):
        zz.append(interpolate_replace_nans(Z[n], kernel))
    fig, axes = plt.subplots(2,1, figsize=(14,20))
    #c = axes[0].contourf(XW, YW, zz_wt, levels=np.linspace(-5, 22, 19), cmap='seismic')
    #c = axes[1].contourf(XM, YM, zz_mut, levels=np.linspace(-5, 22, 19), cmap='seismic')
    
    for n, ax in enumerate(axes.flatten()):
        #m=max_dist+1
        #c= ax.contourf(X, Y, zz[n], levels=np.linspace(0, m, 10), vmin=0, vmax=m,cmap='seismic_r') #not symm
        #f= ax.contour(X, Y, zz[n], levels=np.linspace(0, m, 10), vmin=0, vmax=m,cmap='seismic_r') #not symm
        c= ax.contourf(X, Y, 0.5*(zz[n] + zz[n].T), levels=np.linspace(mn, mx, lev), vmin=mn, vmax=mx, cmap='seismic_r') #symm
        f= ax.contour(X, Y, 0.5*(zz[n] + zz[n].T), levels=np.linspace(mn, mx, lev), vmin=mn, vmax=mx, cmap='seismic_r') #symm
        ax.set_xlim(-10,max(X))
        ax.set_ylim(-10,max(Y))
        ax.set_xticks(np.arange(-10,max(X),5.0))
        ax.set_yticks(np.arange(-10,max(Y),5.0))
        ax.set_xlabel('d1.z ($\AA$)',fontsize=24)
        ax.set_ylabel('d2.z ($\AA$)',fontsize=24)
        ax.set_aspect('equal')
        ax.tick_params(axis='x', labelsize=22)
        ax.tick_params(axis='y', labelsize=22)
        divider=make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad='5%')
        if (n==0):
            cax.set_axis_off()
            ax.set_xlabel('')
        elif (n==1):
            #ax.set_title('del523',fontsize=28, pad=10)
          # cax.set_axis_off()
            #cbar = fig.colorbar(c,cax=cax)
            #cbar.set_label('pka',fontsize=14)
            #cbar.ax.tick_params(labelsize=14)
    
    #fig.subplots_adjust(right=0.8)#  if ((sys=='W') or (prot=='P')):
           cbar = fig.colorbar(c,cax=cax)
           cbar.set_label('distance ($\AA$)',fontsize=24)
          # cbar.set_label('angle ($deg$)',fontsize=24)
           cbar.ax.tick_params(labelsize=22)
    #plt.show()
    #plt.subplots_adjust(hspace=0.2,wspace=0)
    plt.subplots_adjust(hspace=0.2)
    fig.savefig('dist_{nd}_{struc}_{name}.png'.format(nd=nd,struc=struc,name=name))

