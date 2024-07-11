import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import os
from pathlib import Path
import pytraj as pt
import sys

### check and change here ####
##also you may want to change converstion factor between dz data and alldata, colvar and alldata
nframes=int(sys.argv[1])
dmin=int(sys.argv[2])
dmax=int(sys.argv[3])
dx=float(sys.argv[4])
align_mask=':1-510@CA,C,N,O'
system=sys.argv[5]        #as names

if ((system=='wt_Glu0') or (system=='wt_Glu-')):
    strip_mask='!((:1-676))'

elif ((system=='mut_Glu0') or (system=='mut_Glu-')):
    strip_mask='!((:1-675))'

trajspath="/data/mmb/E-Dent/VERONICA/DIMER_LARGER/mw_metad/trajs_mon_wat/phase2"

walkers=np.arange(0,8)
chains=['A','B']
titles=['wild Gext0','wild Gext-','mut Gext0','mut Gext-']
titles2=['WP','WD','MP','MD']
names=['wt_Glu0','wt_Glu-','mut_Glu0', 'mut_Glu-']
###################################################################################################


def select_colvar_bylargest_rbias(colvar_df,system,dvalues,nframes):
    "select nframes of colvar in d1.z and d2.z range with lowest rbias"
    colvar_dz=colvar_df.loc[(colvar_df['system']==system) & (colvar_df['d2.z'] <= dvalues[0]) & 
                    (colvar_df['d2.z'] >= dvalues[1]) & (colvar_df['d1.z'] <= dvalues[2]) &
                    (colvar_df['d1.z'] >= dvalues[3])]
    #print(colvar_dz[:10])
    colvar_largest_rbias=colvar_dz.sort_values(by=['rbias'], ascending=False)[:nframes]
    return colvar_largest_rbias



##select alldata by system, chain, walker, time to get frame of each traj
def select_alldata_bytime(colvar_sel,alldata_df,sys):

    sel_rbias=[]
    selframes=[]
    for ind,row in colvar_sel.iterrows():
        #print(row)
        selframes.append(alldata_df.loc[(alldata_df['system']==row['system']) & (alldata_df['chain']==row['chain'])
        & (alldata_df['walker']==row['walker']) & (alldata_df['time (ps)']==row['time (ps)'])].index.to_list())

        sel_rbias.append(row['rbias'])
    sel_frames=[frame for sublist in selframes for frame in sublist]
        
    return(sel_frames, sel_rbias)




def write_pdb(alldata_df, sel_frames, sel_frames_rbias, trajspath, strip_mask, align_mask):
    for n,fr in enumerate(sel_frames):
        #select traj frame by alldata system, walker,chain and time 
        system=alldata_df['system'].iloc[fr:fr+1].to_string(index=False)
        walker=alldata_df['walker'].iloc[fr:fr+1].to_string(index=False)
        chain=alldata_df['chain'].iloc[fr:fr+1].to_string(index=False)
        traj_frame_from0=int(float(alldata_df['time'].iloc[fr:fr+1].to_string(index=False)))
        topname='{trajspath}/{system}_mon.pdb'.format(trajspath=trajspath,system=system)
        trajname='{trajspath}/{system}_mw{walker}_ch{chain}.xtc'.format(trajspath=trajspath, 
                                                                        system=system, walker=walker,
                                                                        chain=chain)
        print('extracting frame with rbias '+str(sel_frames_rbias[n]))
        current_structure=pt.load(trajname,top=topname,frame_indices=[traj_frame_from0])
        stripped_current=pt.strip(current_structure, strip_mask)
        if (n == 0):
            ref=stripped_current.copy()
            aligned_structure=ref.copy()
        elif (n != 0):
            aligned_structure=pt.align(stripped_current,ref=ref,mask=align_mask, ref_mask=align_mask)
        outnc='{n}.nc'.format(n=n)
        framename='frame{n}'.format(n=n)
        path=Path(framename)
        path.mkdir(parents=True, exist_ok=True)
        os.chdir(path)
        #pt.write_trajectory(outpdb,aligned_structure,format='pdb',overwrite=True)
        pt.write_trajectory(outnc,aligned_structure,overwrite=True)

        #go back (cd ..)
        current_directory = os.getcwd()
        parent_directory = os.path.dirname(current_directory)
        os.chdir(parent_directory)
        ##########################





def write_file(alldata_df, sel_frames, sel_frames_rbias,outfile):
    with open(outfile,'w') as out:
        out.write('PDB index, walker  chain  rbias  d1.z   d2.z\n')
        out.close()
    for n,fr in enumerate(sel_frames):
        system=alldata_df['system'].iloc[fr:fr+1].to_string(index=False)
        walker=alldata_df['walker'].iloc[fr:fr+1].to_string(index=False)
        chain=alldata_df['chain'].iloc[fr:fr+1].to_string(index=False)
        d1z=alldata_df['d1.z'].iloc[fr:fr+1].to_string(index=False)
        d2z=alldata_df['d2.z'].iloc[fr:fr+1].to_string(index=False)
        with open(outfile,'a') as out:
            out.write('{PDB_index},{walker},{chain},{rbias},{d1z},{d2z}\n'.format(PDB_index=n,walker=walker, chain=chain, 
                                                                    rbias=sel_frames_rbias[n],d1z=d1z,d2z=d2z))
    out.close()
    
        
#########################################################################################################################
#                                                    MAIN                                                               #
#########################################################################################################################


#  1  ################################################################
#       LOADING COLVARS DATA TO KNOW RBIAS                           # 
######################################################################

data=[]
path='/data/mmb/E-Dent/MILOSZ/meta/phase2/colvars'
plumed_files=[]
for index,name in enumerate(titles2):
    data.append([])
    for chain in chains:
        for walker in walkers:
            ftemp='{path}/{name}/{w}/COLVAR'.format(path=path,name=name,w=walker)
            if (chain=='A'):
              dtemp=pd.read_csv(ftemp,delimiter=" ",comment='#',skipinitialspace=True,usecols=[0,3,8,25],names=['time','d1.z','d2.z','rbias'])  
            elif (chain=='B'):
              dtemp=pd.read_csv(ftemp,delimiter=" ",comment='#',skipinitialspace=True, usecols=[0,13,18,29],names=['time','d1.z','d2.z','rbias'])
            dtemp['chain']=chain
            dtemp['name']=name
            dtemp['walker']=walker  
            dtemp['time (ps)']=np.round(dtemp['time']).astype(int)
            if (name=='WP'):
                dtemp['system']='wt_Glu0'
            if (name=='WD'):
                dtemp['system']='wt_Glu-'
            if (name=='MP'):
                dtemp['system']='mut_Glu0'
            if (name=='MD'):
                dtemp['system']='mut_Glu-'
            data[index].append(dtemp)
            #plumed_files.append(temp)
colvar_sys=[]
#for each system sys I concatenate the walkers n, poi li appendo. L'ordine è quello di name
for sys,name in enumerate(titles):
    tmp=pd.concat(data[sys][n] for n,m in enumerate(data[sys]))
    colvar_sys.append(tmp)
    
#ora concateno i 4 sistemi    
colvar_allconc=pd.concat((colvar_sys[n] for n,m in enumerate(colvar_sys)),ignore_index=True)
colvar=colvar_allconc.loc[colvar_allconc['time (ps)']%500==0]


#  2  ################################################################
#       LOADING ALLDATA SO THAT I CAN GET                            # 
#       TRAJS INFO (system, walker, chain, frame)                    #
######################################################################


data=[]
for index,name in enumerate(names):
    data.append([])
    for chain in chains:
        for walker in walkers:
            ftemp='/data/mmb/E-Dent/VERONICA/DIMER_LARGER/mw_metad/analysis_phase2/dz/{system}/{system}.mw{w}.{ch}.analysis.dat'.format(system=name,w=walker,ch=chain)
            dtemp=pd.read_csv(ftemp,delimiter=" ",skipinitialspace=True,usecols=[0,1,2],names=['time','d1.z','d2.z'],skiprows=1)
            dtemp['ctrl']=name
            if (name=='wt_Glu0'):
                dtemp['charge']='0'
            if (name=='wt_Glu-'):
                dtemp['charge']='-1'
            if (name=='mut_Glu0'):
                dtemp['charge']='0'
            if (name=='mut_Glu-'):
                dtemp['charge']='-1'
            dtemp['chain']=chain
            dtemp['system']=name
            dtemp['walker']=walker
            dtemp['time (ps)']=(dtemp['time']*500).astype(int)
            data[index].append(dtemp)
data[0][8].head(2)

data_sys=[]
for sys,name in enumerate(names):
    tmp=pd.concat(data[sys][n] for n,m in enumerate(data[sys]))
    data_sys.append(tmp)
alldata=pd.concat((data_sys[n] for n,m in enumerate(data_sys)),ignore_index=True)


#  3  #################################################################
#       -SELECT COLVAR WITH LARGEST RBIAS FOR EACH d1z and d2z range; #
#       -SELECT ALLDATA FRAMES INDEX ( THROUGH COLVAR)                #
#       -RETRIEVE TRAJS INFO FROM ALLDATA AND WRITE PDB & INFO        #
#######################################################################


#if sys dir does not exist, create it!
path=Path(system)
path.mkdir(parents=True, exist_ok=True)
os.chdir(path)

for i in np.linspace(dmin, dmax-dx, int((dmax-dmin)/dx)):
    for j in np.linspace(dmin, dmax-dx, int((dmax-dmin)/dx)): 
        print(i)
        print(j)
        #it can be chainA or chainB
        colvar_largest_rbias=select_colvar_bylargest_rbias(colvar,system,[j+dx,j,i+dx,i],nframes)
        sel_frames, sel_frames_rbias=select_alldata_bytime(colvar_largest_rbias,alldata,system)
        print(sel_frames)
        outdir='frames'+'_'+str(i)+'_'+str(j) #i=dz1, j=dz2
        outfile='frames_{dz1}_{dz2}.dat'.format(dz1=i,dz2=j)
        #if bin dir doesn't exist, create it!
        path=Path(outdir)
        path.mkdir(parents=True, exist_ok=True)
        os.chdir(path)
        
        write_file(alldata,sel_frames,sel_frames_rbias,outfile)
        write_pdb(alldata,sel_frames,sel_frames_rbias,trajspath,strip_mask,align_mask)
        os.chdir("..")

