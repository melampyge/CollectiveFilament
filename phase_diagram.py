#!/usr/bin/python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
import seaborn as sns
import os.path
import glob
import pandas as pd
from string import atof
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import colorsys 
import seaborn as sns
sns.set()
sns.set(style="white",context='notebook',
        font_scale=1.4,font="Open Sans",
        rc={'mathtext.default': 'regular','font.size': 20, 
            'font.family': 'sans',"figure.dpi":120,"axes.formatter.limits":(-2,2)})


def plot_slices(Data,sliceby,namex,namey,save=False,scale=.15,savesuffix='',binning=False):
    labels={'Pe':r"$Pe = f_pL^2/k_BT$",
           'xi_L':r"$\xi_p/L = \kappa/k_BTL$",
           'density':'Density'}
    
    if isinstance(binning,int) and binning:
        bins = np.linspace(Data[sliceby].min(), Data[sliceby].max(), binning)
        Data=Data.copy()
        Data['digi']=Data[sliceby].apply(lambda x:np.digitize(x,bins))
        group = Data.groupby('digi')
    else:
        group=Data.groupby(sliceby)

    for i,grp in group:
        
        grp=grp.copy()
        fwidth=10.
        fig,ax=plt.subplots(1,1,figsize=(fwidth,fwidth*3/4.))
        ax.set_title(r'{1:s}: {0:.2f}'.format(grp[sliceby].mean(),labels[sliceby]))

        namex2x=dict([(j,i) for i,j in enumerate(sorted(grp[namex].unique()))])
        namey2y=dict([(j,i) for i,j in enumerate(sorted(grp[namey].unique()))])
        grp['x']=grp[namex].apply(lambda x:namex2x[x])
        grp['y']=grp[namey].apply(lambda x:namey2y[x])


        ax.scatter(grp.x,grp.y)
        ax.set_xlabel(labels[namex])
        ax.set_ylabel(labels[namey])

        show_screenshot=True
        if show_screenshot:
            for fname,x,y in zip(grp.folder,grp.x,grp.y):

                if os.path.exists(fname):
                    arr_hand = read_png(fname)
                
                    zoom=scale
                    if arr_hand.shape[0]==1024:
                        zoom/=2.
                    imagebox = OffsetImage(arr_hand, zoom=zoom)
    
                    xy = [x,y]               # coordinates to position this image
    
                    ab = AnnotationBbox(imagebox, xy,
                        xybox=(0., -0.),
                        xycoords='data',
                        boxcoords="offset points",frameon=1,pad=.1)                                  
                    ax.add_artist(ab)
    
            ax.xaxis.set_ticks(sorted(grp.x.unique()))
            ax.yaxis.set_ticks(sorted(grp.y.unique()))
    
            xticks=sorted(grp[namex].unique().tolist())
            yticks=sorted(grp[namey].unique().tolist())
            xlabels=["{0:.1f}".format(j) for j in xticks]
            ylabels=["{0:.2f}".format(j) for j in yticks]
            ax.xaxis.set_ticklabels(xlabels)
            ax.yaxis.set_ticklabels(ylabels)
    
            ax.grid(1,color='#cccccc',linestyle='--')
            if save:
                fig.savefig("many_polymers_5/plots/set_{1:s}{0:.2f}.pdf".format(grp[sliceby].mean(),sliceby+savesuffix),bbox_inches='tight',dpi=300)    


## Index the data
density = [0.08]
kappa = [2.5, 5.0, 25.0, 62.5, 125.0, 200.0]
fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]

folders = []
for d in density:
    for k in kappa:
        for f in fp:
            folders.append( 'many_polymers_5/density_' + str(d) + '/kappa_' + str(k) 
            + '/fp_' + str(f) )


## Manually read from the input files
kbT = 1.
L = 25.0

print np.size(folders)
files = []
tmp = {}
for folder in folders:
    
    tmp = {}
    dname = folder.split('/')[1]
    ddict = dict(zip(dname.split('_')[::2],map(atof,dname.split('_')[1::2])))
    
    kname = folder.split('/')[2]
    kdict = dict(zip(kname.split('_')[::2],map(atof,kname.split('_')[1::2])))
    
    fname = folder.split('/')[3]
    fdict = dict(zip(fname.split('_')[::2],map(atof,fname.split('_')[1::2])))

    tmp.update(ddict)
    tmp.update(kdict)
    tmp.update(fdict)

    peclet = tmp['fp']*L**2/kbT
    persistence_over_L = tmp['kappa']/(L*kbT)    
    
    tmp.update({'folder':folder+'/fin.png', 'Pe':peclet, 'xi_L':persistence_over_L})
    files.append(tmp)

Data = pd.DataFrame(files,columns=files[0].keys())
sliceby = 'Pe'
bins = np.linspace(Data[sliceby].min(), Data[sliceby].max(), 5)

print files
print files[0]

plot1={'sliceby':'density',
'namex':'Pe',
'namey':'xi_L'}


plot_slices(Data,save=True,**plot1)


