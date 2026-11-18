#!/usr/bin/env python3
"""Current-cohort candidate references; reproduce with python compare_models.py."""
from pathlib import Path
import hashlib, json, time
import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=ROOT/'outputs/yaw/models'
SOURCE=ROOT/'outputs/yaw'
MODELS=['current','angular','target_cv','spatial']
LABELS=['Current target aligned','Constant angular lead','Constant-velocity extrapolation','Constant spatial lead']
COLORS=['#ABB0B5','#E69F00','#7E2F8E','#009E73']
GRIDS=[np.array([0.]),np.deg2rad(np.arange(0,180.001,.5)),np.arange(0,.40001,.005),np.arange(0,5.00001,.025)]
UNITS=['none','radians','seconds','metres']
DEG=180/np.pi
CONTEXTS=['frontal pursuit','lateral pursuit','close pass-by pursuit']
def wrap(x): return (x+np.pi)%(2*np.pi)-np.pi
def yaw(p):
    x,y,z=p[...,0],p[...,1],p[...,2]
    return np.arctan2(-np.sin(z)*np.cos(y)-np.cos(z)*np.sin(x)*np.sin(y),np.cos(z)*np.cos(y)-np.sin(z)*np.sin(x)*np.sin(y))
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def velocity(x,frames,fs):
    n=max(3,int(np.floor(.1*fs+.5))+1)
    t=np.arange(1-n,1)/fs
    coef=np.linalg.pinv(np.c_[np.ones(n),t,t*t])[1]
    out=np.full_like(x,np.nan)
    for end in range(n-1,len(frames)):
        start=end-n+1
        if np.all(np.diff(frames[start:end+1])==1):
            out[:,end]=np.einsum('k,dkc->dc',coef,x[:,start:end+1])
    return out

def states(row,accepted):
    name=row.dataset; folder=Path(row.dataDir)
    paths=[folder/'MonteCarloResults'/f'mcMain_{name}.mat',folder/'storedOptimRaw_trajectory.mat',folder/'3DLurePos.mat',folder/'lureAnalysis_results_trajectory.mat']
    r=loadmat(paths[3],simplify_cells=True)['runResults']; fs=float(r['fs'])
    af=accepted[accepted.dataset==name].absFrame.to_numpy(dtype=int)
    with h5py.File(paths[0]) as h:
        allframes=np.asarray(h['absFrames']).ravel().astype(int)
        ix=np.array([np.flatnonzero(allframes==f).item() for f in af])
        ok=np.flatnonzero(np.asarray(h['sampleSucceeded']).ravel())
        assert len(ok)==500
        nominal=np.asarray(h['nominalPose']).T[ix]
        lure=np.asarray(h['nominalLureXYZ']).T[ix]
        pose=np.concatenate([nominal[None],np.asarray(h['poseAll']).transpose(2,1,0)[ok][:,ix]])
        target=np.concatenate([lure[None],np.asarray(h['lureXYZAll']).transpose(2,1,0)[ok][:,ix]])
    # Verify current sources, not only historical file names.
    p=loadmat(paths[1],simplify_cells=True)
    pi=np.array([np.flatnonzero(np.asarray(p['absFrames']).ravel()==f).item() for f in af])
    assert np.allclose(nominal,np.asarray(p['xBarRawStore'])[pi],atol=1e-9,rtol=0),name
    l=loadmat(paths[2],simplify_cells=True)['Lure3D']
    li=np.array([np.flatnonzero(np.asarray(l['Index']).ravel()==f).item() for f in af])
    marks=np.asarray(l['Mark'])
    if marks.dtype==object: marks=np.vstack([np.asarray(z).ravel()[:3] for z in marks.ravel()])
    assert np.allclose(lure,marks[li,:3],atol=1e-9,rtol=0),name
    ri=np.array([np.flatnonzero(np.asarray(r['commonFrames']).ravel()==f).item() for f in af])
    head=yaw(pose[:,:,:3]); torso=yaw(pose[:,:,6:9])
    rel=target[:,:,:2]-pose[:,:,3:5]; rng=np.linalg.norm(rel,axis=-1)
    lam=np.arctan2(rel[:,:,1],rel[:,:,0])
    for key,z in [('theta_h_yaw',head),('theta_t_yaw',torso),('lambda_h_yaw',lam)]:
        assert np.max(np.abs(wrap(np.asarray(r[key])[ri]-z[0])))<1e-8,(name,key)
    assert np.asarray(r['yawAnalysisMask'],bool)[ri].all(),name
    v=velocity(target[:,:,:2],af,fs); speed=np.linalg.norm(v,axis=-1)
    valid=np.all(np.isfinite(head)&np.isfinite(torso)&np.isfinite(speed)&(speed>1e-6)&(rng>1e-4),axis=0)
    assert valid.sum()>2,name
    # QC diagnostics use exactly the shared model evaluation frames.
    bearing_dev=wrap(lam[1:]-lam[0])*DEG
    bearing_width=np.percentile(bearing_dev,97.5,axis=0)-np.percentile(bearing_dev,2.5,axis=0)
    td=wrap(torso[0]-lam[0]); hd=wrap(head[0]-lam[0])
    between=(np.abs(hd)<=np.abs(td)+1e-12)&(hd*td>=-1e-12)
    q=dict(dataset=name,animal=row.animal,context=row.context,acceptedFrames=len(af),scoredFrames=int(valid.sum()),
           fsHz=fs,rangeMedianM=float(np.median(rng[0,valid])),rangeMinM=float(np.min(rng[0,valid])),
           speedMedianMps=float(np.median(speed[0,valid])),bearing95WidthMedianDeg=float(np.median(bearing_width[valid])),
           forecast50msToRangeMedian=float(np.median(.05*speed[0,valid]/rng[0,valid])),
           headBetweenTargetTorsoPercent=float(100*between[valid].mean()))
    return dict(head=head[:,valid],torso=torso[:,valid],lam=lam[:,valid],rel=rel[:,valid],v=v[:,valid],
                speed=speed[:,valid],frames=af[valid],qc=q,paths=paths)

def score_curves(s):
    losses=[]; biases=[]; singular=[]
    for model,grid in zip(MODELS,GRIDS):
        loss=np.empty((501,len(grid))); bias=np.empty_like(loss); near=0
        for start in range(0,501,20):
            sl=slice(start,min(start+20,501));lam=s['lam'][sl];rel=s['rel'][sl];v=s['v'][sl]
            if model=='current': aim=lam[:,:,None]
            elif model=='angular':
                cross=rel[:,:,0]*v[:,:,1]-rel[:,:,1]*v[:,:,0]
                sign=np.sign(cross);sign[np.abs(cross)<1e-12*np.linalg.norm(rel,axis=-1)*s['speed'][sl]]=0
                aim=lam[:,:,None]+sign[:,:,None]*grid
            elif model in ['target_cv','spatial']:
                vector=v if model=='target_cv' else v/s['speed'][sl,:,None]
                pos=rel[:,:,None,:]+vector[:,:,None,:]*grid[None,None,:,None]
                near+=int((np.linalg.norm(pos,axis=-1)<1e-6).sum())
                aim=np.arctan2(pos[:,:,:,1],pos[:,:,:,0])
            else: raise ValueError(model)
            err=wrap(s['head'][sl,:,None]-aim)
            loss[sl]=np.mean(err**2,axis=1);bias[sl]=np.mean(err,axis=1)**2
        assert near==0,(model,near)
        losses.append(loss);biases.append(bias);singular.append(near)
    return losses,biases

def weights(cohort,ix,scheme):
    if scheme=='LORO': return np.ones(len(ix))/len(ix)
    animals=cohort.iloc[ix].animal.to_numpy();names,counts=np.unique(animals,return_counts=True)
    return np.array([1/(len(names)*counts[np.flatnonzero(names==a)[0]]) for a in animals])

def compare(cohort,losses,biases,ix,scheme,transfer=False):
    pred=np.full((501,len(MODELS),len(ix)),np.nan); bias=np.full_like(pred,np.nan); pars=np.full_like(pred,np.nan)
    for k,test in enumerate(ix):
        train=[i for i in (list(np.flatnonzero(cohort.context=='lateral pursuit')) if transfer else ix)
               if i!=test and (scheme=='LORO' or cohort.iloc[i].animal!=cohort.iloc[test].animal)]
        assert train
        w=weights(cohort,train,scheme)
        for m in range(len(MODELS)):
            objective=np.einsum('r,rdp->dp',w,losses[m][train])
            best=objective.argmin(axis=1)
            pred[:,m,k]=losses[m][test,np.arange(501),best]
            bias[:,m,k]=biases[m][test,np.arange(501),best]
            pars[:,m,k]=GRIDS[m][best]
    return pred,bias,pars

def interval(x):
    q=np.percentile(x[1:],[50,2.5,97.5]);return dict(nominal=float(x[0]),median=float(q[0]),lower=float(q[1]),upper=float(q[2]))

def plot_bars(rows,groups,path,aggregate_label):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(max(11.5,len(groups)*1.8),6.5))
    fig.subplots_adjust(left=.075,right=.99,bottom=.16,top=.84)
    lateral_only = path.stem == 'candidate_model_comparison_lateral5'
    omit_torso = lateral_only or path.stem == 'candidate_model_comparison_all18'
    models = MODELS[:4] if omit_torso else MODELS
    colors = COLORS[:4] if omit_torso else COLORS
    labels = ['Current target aligned'] + LABELS[1:4] if omit_torso else LABELS
    width=.18 if omit_torso else .15
    centre=(len(models)-1)/2
    x=np.arange(len(groups),dtype=float); x[-1]+=.45; tops=[]
    for m,(model,color) in enumerate(zip(models,colors)):
        z=rows.set_index(['group','model']).loc[[(g,model) for g in groups]]
        mid=z['median'].to_numpy();lo=z.lower.to_numpy();hi=z.upper.to_numpy();pos=x+(m-centre)*width
        pos[-1]=x[-1]+(m-centre)*.22
        ax.bar(pos,mid,width=width,color=color,zorder=2)
        ax.errorbar(pos,mid,yerr=np.array([mid-lo,hi-mid]),fmt='none',ecolor='black',elinewidth=1.2,capsize=3,zorder=3)
        tops.extend(hi)
        for j in (range(len(groups)) if lateral_only else [len(groups)-1]):
            ax.annotate(f'{mid[j]:.1f}°',(pos[j],hi[j]),xytext=(0,6+12*(m%2)),textcoords='offset points',ha='center',fontsize=8.5)
    ax.set_ylabel('Head-yaw prediction error (RMSE, °)')
    tick_labels = [f'Run {j+1}' for j in range(len(groups)-1)] + [aggregate_label] if lateral_only else [g.replace(' pursuit','\npursuit') if g!=groups[-1] else aggregate_label for g in groups]
    ax.set_xticks(x,tick_labels)
    ax.set_ylim(0,max(tops)*1.18+2)
    ax.axvline((x[-2]+x[-1])/2,color='.65',ls='--',lw=1)
    ax.spines[['top','right']].set_visible(False)
    legend=[Patch(color=c,label=l) for c,l in zip(colors,labels)]
    dummy=ax.errorbar([np.nan],[np.nan],yerr=[1],fmt='none',ecolor='black',capsize=3,label='MC 95% reconstruction interval')
    fig.legend(handles=legend+[dummy],loc='upper center',bbox_to_anchor=(.53,.995),ncol=3,frameon=True,fancybox=False,edgecolor='.25',fontsize=10)
    for ext in ['png','pdf','svg']:fig.savefig(path.with_suffix('.'+ext),dpi=400,facecolor='white')
    plt.close(fig)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cohort=pd.read_csv(SOURCE/'yaw_absolute_error_cohort.csv')
    cohort=cohort[cohort.context=='lateral pursuit'].reset_index(drop=True)
    accepted=pd.read_csv(SOURCE/'yaw_absolute_error_by_frame.csv')
    accepted=accepted[accepted.dataset.isin(cohort.dataset)]
    assert len(cohort)==5 and len(accepted)==510
    all_loss=[[] for _ in MODELS]; all_bias=[[] for _ in MODELS]
    qc=[]; scored=[]; manifest=[]; started=time.monotonic()
    for _,row in cohort.iterrows():
        s=states(row,accepted); lc,bc=score_curves(s)
        for m in range(len(MODELS)): all_loss[m].append(lc[m]); all_bias[m].append(bc[m])
        qc.append(s['qc'])
        scored.extend(dict(dataset=row.dataset,absFrame=int(f)) for f in s['frames'])
        manifest.extend(dict(path=str(p),sha256=digest(p)) for p in s['paths'])
        print(row.dataset,s['qc']['scoredFrames'],'evaluation frames',flush=True)
    losses=[np.stack(z) for z in all_loss]; biases=[np.stack(z) for z in all_bias]
    pd.DataFrame(qc).to_csv(OUT/'geometry_diagnostics.csv',index=False)
    pd.DataFrame(scored).to_csv(OUT/'scored_frames.csv',index=False)
    summary=[]; by_run=[]; parameters=[]; draws={}; ix=list(range(len(cohort)))
    for scheme in ['LORO','LOAO']:
        val,bias,pars=compare(cohort,losses,biases,ix,scheme)
        agg=np.sqrt(np.einsum('r,dmr->dm',weights(cohort,ix,scheme),val))*DEG
        draws[scheme]=val; draws[scheme+'_parameters']=pars
        for m,model in enumerate(MODELS):
            summary.append(dict(cohort='lateral5',scheme=scheme,model=model,frames=len(scored),**interval(agg[:,m])))
            for j,row in cohort.iterrows():
                by_run.append(dict(cohort='lateral5',scheme=scheme,model=model,dataset=row.dataset,
                                   animal=row.animal,**interval(np.sqrt(val[:,m,j])*DEG)))
                parameters.append(dict(scheme=scheme,model=model,dataset=row.dataset,unit=UNITS[m],**interval(pars[:,m,j])))
    for name,table in [('summary',summary),('by_run',by_run),('parameters',parameters)]:
        pd.DataFrame(table).to_csv(OUT/f'{name}.csv',index=False)
    np.savez_compressed(OUT/'heldout_draws.npz',**draws)
    rows=[]; groups=[]
    for j,row in cohort.iterrows():
        label=f'L{j+1}'; groups.append(label)
        rows.extend(dict(z,group=label) for z in by_run if z['scheme']=='LORO' and z['dataset']==row.dataset)
    groups.append('Summary')
    rows.extend(dict(z,group='Summary') for z in summary if z['scheme']=='LORO')
    plot_bars(pd.DataFrame(rows),groups,OUT/'candidate_model_comparison_lateral5','Equal-run\nsummary')
    (OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'verification.json').write_text(json.dumps(dict(acceptedFrames=510,scoredFrames=len(scored),
        runs=5,draws=500,seconds=time.monotonic()-started,validation=['LORO','LOAO']),indent=2)+'\n')
    print(pd.DataFrame(summary).to_string(index=False))

if __name__=='__main__': main()
