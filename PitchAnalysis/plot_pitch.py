"""Render the refreshed MATLAB statistics without headless MATLAB font overlap."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
B=Path(__file__).resolve().parents[1]/'outputs'/'pitch'; T=B/'tightness_run_equal'; M=B/'mechanism'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'legend.framealpha':1,'legend.fancybox':False,'legend.edgecolor':'.4','pdf.fonttype':42})
H=(.38,.40,.88); O=(.95,.52,.28); G=(.20,.62,.35)
def axes(size=(7.6,5.6)):
    f,a=plt.subplots(figsize=size,layout='constrained'); a.grid(alpha=.18); a.set_axisbelow(True); return f,a
def save(f,d,n):
    for ext in ('png','pdf'): f.savefig(d/f'{n}.{ext}',dpi=300,facecolor='white')
    plt.close(f)
frames=pd.read_csv(T/'pitch_frame_level_errors.csv'); assert len(frames)==1590 and frames.datasetKey.nunique()==16
r=pd.read_csv(T/'pitch_absolute_error_ecdf_by_run.csv'); e=pd.read_csv(T/'pitch_absolute_error_ecdf.csv')
f,a=axes()
handles=[]
for col,mean,se,color,name in [('headCdf','headMeanCdf','headSE',H,'Head'),('torsoCdf','torsoMeanCdf','torsoSE',O,'Torso')]:
    curves=[]
    for key,v in r.groupby('datasetKey',sort=False):
        x=v.absoluteErrorDeg.to_numpy(); y=v[col].to_numpy(); curves.append(y)
        values=frames[frames.datasetKey==key]['absAlphaDeg' if name=='Head' else 'absBetaDeg'].to_numpy()
        np.testing.assert_allclose(y,(values[:,None]<=x).mean(axis=0),atol=1e-12)
        a.plot(x,y,color=color,alpha=.22,lw=.7)
    np.testing.assert_allclose(e[mean],np.mean(curves,axis=0),atol=1e-12)
    np.testing.assert_allclose(e[se],np.std(curves,axis=0,ddof=1)/np.sqrt(16),atol=1e-12)
    band=a.fill_between(e.absoluteErrorDeg,np.maximum(0,e[mean]-e[se]),np.minimum(1,e[mean]+e[se]),color=color,alpha=.13,label=f'{name}: ±1 SE')
    line,=a.plot(e.absoluteErrorDeg,e[mean],color=color,lw=2.5,label=f'{name}–lure pitch error')
    handles.extend([line,band])
a.axvline(10,color='.3',ls='--',lw=1);a.text(9.7,.025,'10°',rotation=90,ha='right')
a.set(xlim=(0,e.absoluteErrorDeg.max()),ylim=(0,1),xlabel='Absolute pitch error (deg)',ylabel='Cumulative proportion of frames');a.legend(handles=handles,loc='lower right',fontsize=10)
save(f,T,'pitch_absolute_error_ecdf')
c=pd.read_csv(T/'pitch_centered_error_distribution.csv'); x=c.centeredErrorDeg.to_numpy();bw=np.diff(x)[0];edges=np.r_[x-bw/2,x[-1]+bw/2]
f,a=axes(); hs=[]
for col,mean,se,color,name in [('alphaDeg','headMeanProp','headSE',H,'Head'),('betaDeg','torsoMeanProp','torsoSE',O,'Torso')]:
    curves=[]
    for _,v in frames.groupby('datasetKey',sort=False):
        z=v[col].to_numpy();hist=np.histogram(z-z.mean(),edges)[0];p=hist/hist.sum();curves.append(p);a.plot(x,p,color=color,alpha=.2,lw=.6)
    np.testing.assert_allclose(c[mean],np.mean(curves,axis=0),atol=1e-12)
    np.testing.assert_allclose(c[se],np.std(curves,axis=0,ddof=1)/np.sqrt(16),atol=1e-12)
    y=c[mean].to_numpy();band=a.fill_between(x,np.maximum(0,y-c[se]),np.minimum(1,y+c[se]),color=color,alpha=.13,label=f'{name}: ±1 SE across runs')
    l,=a.plot(x,y,color=color,lw=2.5,label=f'{name}–lure pitch error');hs.extend([l,band])
    mask=np.isin(x,np.arange(-15,16,5));a.errorbar(x[mask],y[mask],yerr=c[se][mask],fmt='o',ms=5,color=color,mec='white',capsize=3)
    peak=np.argmax(y);half=y.max()/2
    left=np.interp(half,y[:peak+1],x[:peak+1]);right=np.interp(half,y[peak:][::-1],x[peak:][::-1]);a.plot([left,right],[half,half],color=color,lw=1.2)
a.axvline(0,color='.5',ls=':',lw=1);a.set(xlim=(x.min()-bw,x.max()+bw),ylim=(0,None),xticks=x,xlabel='Run-centred pitch error (deg)',ylabel='Proportion of frames');a.legend(handles=hs,loc='upper right',fontsize=10)
save(f,T,'pitch_centered_error_distribution')
r=pd.read_csv(M/'pitch_control_model_by_run.csv');r=r[np.isclose(r.requestedCorrectionDelaySec,.05)];assert len(r)==16
p=pd.read_csv(M/'pitch_continuous_disturbance_rejection.csv');assert p.datasetKey.nunique()==16
lim=max(10,np.ceil(np.percentile(np.abs(p[['opposingTorsoChangeDeg','qChangeDeg']].to_numpy()),99)/5)*5);edges=np.linspace(-lim,lim,51);hist=[]
for _,v in p.groupby('datasetKey'):
    z=np.histogram2d(v.opposingTorsoChangeDeg,v.qChangeDeg,bins=(edges,edges))[0];hist.append(z/z.sum())
density=np.mean(hist,axis=0);assert np.isclose(density.sum(),1)
f,a=axes((6.5,6));a.grid(False);im=a.imshow(density.T,origin='lower',extent=(-lim,lim,-lim,lim),cmap=LinearSegmentedColormap.from_list('green',['white',G]),aspect='equal');cb=f.colorbar(im,ax=a,shrink=.8);cb.set_label('Run-equal probability per bin')
g=r.counterRotationDisturbanceGain; se=g.std(ddof=1)/np.sqrt(len(g));a.plot([-lim,lim],[-lim,lim],'--',color='.3',lw=1.3,label='Complete cancellation ($g_D=1$)');a.plot([-lim,lim],np.array([-lim,lim])*g.mean(),color=G,lw=2,label=f'Run-equal gain: $g_D={g.mean():.2f}\\pm{se:.2f}$')
a.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel=r'Opposing torso pitch change, $-\Delta\theta_T$ (deg)',ylabel=r'Head-on-torso pitch change, $\Delta q$ (deg)');a.legend(loc='upper left',fontsize=9)
save(f,M,'pitch_continuous_disturbance_rejection_run_equal')
n_animals=r.animal.nunique(); pop_x=n_animals+1.15
f,a=axes((8,5.6));colors=['#0072bd','#d95319','#edb120','#7e2f8e','#77ac30','#4dbeee'];means=[]
for i,(animal,v) in enumerate(r.groupby('animal',sort=False),1):
    y=v.torsoTransmissionGain.to_numpy();mu=y.mean();means.append(mu);xs=i+np.linspace(-.17,.17,len(y));color=colors[i-1]
    for xx,yy in zip(xs,y):a.plot([xx,i],[yy,mu],color=color,alpha=.32,lw=.8)
    a.scatter(xs,y,s=30,color=color,edgecolor='white',zorder=3);a.scatter([i],[mu],s=120,color=color,edgecolor='.1',zorder=4)
mu=np.mean(means);se=np.std(means,ddof=1)/np.sqrt(len(means));a.errorbar(pop_x,mu,yerr=se,fmt='D',mfc='white',color='.1',ms=8,capsize=5,lw=1.8,zorder=4)
a.axhline(0,color='.2',ls='--');a.axhline(1,color='.45',ls='--');a.text(.55,.035,'$g_T=0$: complete rejection',fontsize=10);a.text(.55,1.035,'$g_T=1$: complete transmission',fontsize=10)
a.set(xlim=(.45,pop_x+.95),ylim=(-.55,1.15),xticks=list(range(1,n_animals+1))+[pop_x],xticklabels=[f'Cheetah {i}' for i in range(1,n_animals+1)]+['Pop.'],ylabel='$g_T$');a.tick_params(axis='x',labelsize=9);a.grid(axis='x',visible=False)
a.legend(handles=[Line2D([],[],color='.6',lw=.8,label='Run-to-animal mean'),Line2D([],[],marker='o',ls='',color='.4',ms=5,label='Individual run'),Line2D([],[],marker='o',ls='',color='.4',ms=10,label='Animal mean'),Line2D([],[],marker='D',ls='',color='.1',mfc='white',label='Population mean ±1 SE')],loc='upper right',bbox_to_anchor=(1,.83),fontsize=9)
save(f,M,'pitch_population_control_gains_run_level')
t=pd.read_csv(M/'pitch_representative_control_timeseries.csv');f,a=axes((11.2,5))
for col,color,style,label in [('thetaTdeg',(.88,.30,.07),'-',r'Torso pitch, $\theta_T$'),('qDeg',G,'-',r'Head-on-torso pitch, $q$'),('thetaHdeg',(.05,.43,.72),'-',r'Head pitch, $\theta_H$'),('lambdaDeg',(.12,.12,.12),'--',r'Target elevation, $\lambda$'),('alphaDeg',(.55,.25,.65),':',r'Head–target error, $\alpha$')]:a.plot(t.timeSec,t[col],color=color,ls=style,lw=1.8,label=label)
a.set(xlim=(.2,t.timeSec.max()),xlabel='Time (s)',ylabel='Pitch angle (deg)')
a.legend(loc='lower center',bbox_to_anchor=(.5,1.02),ncol=3,fontsize=10)
f.suptitle('Representative pursuit sequence: phantomTopRun1_20170903',fontsize=12)
save(f,M,'pitch_representative_control_timeseries')
print('Verified: 16 runs, 1590 frames; ECDFs and centred distributions reproduce MATLAB means and SEs; density sums to one.')
