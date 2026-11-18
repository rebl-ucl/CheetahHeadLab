#!/usr/bin/env python3
"""Uncentred signed gamma/delta distributions, equal weight per lateral pursuit."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from plotYawAbsoluteErrorEcdf import HEAD,TORSO

OUT=Path(__file__).resolve().parents[1]/'outputs'/'yaw'
STEM='yaw_signed_error_distribution_lateral_only'

def main(all_runs=False, include_close=False):
    assert not (all_runs and include_close)
    stem = "yaw_signed_error_distribution_all18" if all_runs else ("yaw_signed_error_distribution_yaw_only" if include_close else STEM)
    source=OUT/'yaw_absolute_error_by_frame.csv'
    frames=pd.read_csv(source)
    if not all_runs:
        frames=frames.loc[frames.context.isin(['lateral pursuit','close pass-by pursuit'] if include_close else ['lateral pursuit'])].copy()
    n_runs=18 if all_runs else (7 if include_close else 5)
    n_frames=1708 if all_runs else (628 if include_close else 510)
    n_animals=6 if all_runs else (4 if include_close else 3)
    assert len(frames)==n_frames and frames.dataset.nunique()==n_runs and frames.animal.nunique()==n_animals
    if all_runs:
        assert frames.groupby('context').size().to_dict()=={'frontal pursuit':1080,'lateral pursuit':510,'close pass-by pursuit':118}
    if include_close:
        assert frames.groupby('context').size().to_dict()=={'lateral pursuit':510,'close pass-by pursuit':118}
    columns=['headMedianYawDeg','torsoMedianYawDeg']
    values=frames[columns].to_numpy()
    assert np.isfinite(values).all() and ((values>=-180)&(values<180)).all()
    bw=5.
    lo=min(0,np.floor(values.min()/bw)*bw)-bw
    hi=max(0,np.ceil(values.max()/bw)*bw)+bw
    x=np.arange(lo,hi+bw/2,bw);edges=np.r_[x-bw/2,x[-1]+bw/2]
    table=pd.DataFrame({'signedYawErrorDeg':x});byrun=[];curves={};stats={}
    for segment in ['head','torso']:
        arr=[];means=[];variances=[]
        for dataset,g in frames.groupby('dataset',sort=False):
            raw=g[segment+'MedianYawDeg'].to_numpy()
            counts=np.histogram(raw,bins=edges)[0]
            assert counts.sum()==len(g)
            prop=counts/len(g);arr.append(prop);means.append(raw.mean());variances.append(raw.var())
            byrun.extend(dict(dataset=dataset,segment=segment,signedYawErrorDeg=a,proportion=b) for a,b in zip(x,prop))
        arr=np.asarray(arr);curves[segment]=arr
        mean=arr.mean(axis=0);se=arr.std(axis=0,ddof=1)/np.sqrt(n_runs)
        assert np.isclose(mean.sum(),1)
        table[segment+'MeanProp']=mean;table[segment+'SE']=se
        stats[segment]={'equalRunMeanDeg':float(np.mean(means)),'equalRunMixtureSdDeg':float(np.sqrt(np.mean(variances)+np.var(means)))}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(7.6,5.6));fig.subplots_adjust(left=.12,right=.98,bottom=.14,top=.97)
    ymax=0
    for segment,color in [('torso',TORSO),('head',HEAD)]:
        mean=table[segment+'MeanProp'].to_numpy();se=table[segment+'SE'].to_numpy()
        ax.plot(x,curves[segment].T,color=color,alpha=.20,lw=.6)
        ax.fill_between(x,np.maximum(0,mean-se),np.minimum(1,mean+se),color=color,alpha=.13,lw=0)
        ax.plot(x,mean,color=color,lw=2.5)
        mask=mean>0
        ax.errorbar(x[mask],mean[mask],yerr=se[mask],fmt='o',color=color,mfc=color,mec='white',ms=4,lw=1,capsize=2.5)
        ymax=max(ymax,curves[segment].max(),(mean+se).max())
    ax.axvline(0,color='.5',ls=':',lw=.8)
    ax.set(xlabel='Signed yaw error (deg)',ylabel='Proportion of frames',xlim=(lo,hi),ylim=(-.02,ymax*1.08))
    ax.grid(alpha=.35,lw=.6);ax.tick_params(direction='in')
    ax.legend(handles=[Line2D([],[],color=HEAD,lw=2.5,label='Head–lure yaw error'),Patch(facecolor=HEAD,alpha=.13,edgecolor='none',label='Head: ±1 SE across runs'),Line2D([],[],color=TORSO,lw=2.5,label='Torso–lure yaw error'),Patch(facecolor=TORSO,alpha=.13,edgecolor='none',label='Torso: ±1 SE across runs')],loc='upper right',frameon=True,fancybox=False,edgecolor='.25',fontsize=9)
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{stem}.{ext}',dpi=600,facecolor='white')
    plt.close(fig)
    table.to_csv(OUT/f'{stem}.csv',index=False)
    pd.DataFrame(byrun).to_csv(OUT/f'{stem}_by_run.csv',index=False)
    frames.to_csv(OUT/f'{stem}_frames.csv',index=False)
    caption=('Uncentred signed head–lure (gamma, blue) and torso–lure (delta, orange) yaw-error distributions across five lateral pursuits (510 frames, three cheetahs). Original wrapped posterior-median signed errors are retained, without subtracting run means, taking absolute values or reorienting signs by target travel. Bins are 5 degrees wide. Thick lines show equally weighted means of per-run proportions; thin lines show individual runs; shading and error bars show +/-1 SE across the five runs. All frames are included. Zero denotes current-target alignment; the sign is the original angular-coordinate sign, not the travel-directed ahead/lagged classification.\n')
    if all_runs:
        caption=caption.replace('five lateral pursuits (510 frames, three cheetahs)', '18 pursuits (1,708 frames, six cheetahs): 11 frontal (1,080 frames), five lateral (510 frames) and two close pass-by pursuits (118 frames)').replace('across the five runs', 'across the 18 runs')
        caption += 'Each run has equal weight; contexts therefore contribute 11/18, 5/18 and 2/18 of the aggregate. Differences in width include both within-run variation and between-run/context offsets.\n'
    if include_close:
        caption=caption.replace('five lateral pursuits (510 frames, three cheetahs)', 'seven yaw-tracking pursuits (628 frames, four cheetahs): five lateral pursuits (510 frames) and two close pass-by pursuits (118 frames)').replace('across the five runs','across the seven runs')
    (OUT/f'{stem}_caption.txt').write_text(caption)
    report={'nRuns':n_runs,'nAnimals':n_animals,'nFrames':n_frames,'binWidthDeg':bw,'centred':False,'weighting':'equal run','statistics':stats,'source':str(source),'sourceSha256':hashlib.sha256(source.read_bytes()).hexdigest(),'validation':'All original signed errors retained unchanged; all frames binned; per-run and aggregate distributions sum to one.'}
    (OUT/f'{stem}_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-runs',action='store_true',help='Use all 18 pursuits instead of five lateral pursuits.')
    parser.add_argument('--include-close',action='store_true',help='Include two close pass-by pursuits with the five lateral pursuits.')
    args=parser.parse_args()
    main(args.all_runs,args.include_close)
