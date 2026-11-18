#!/usr/bin/env python3
"""Reproduce the seven-yaw-run subset using the accepted paired-frame tables."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import plotYawAbsoluteErrorEcdf as renderer

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'outputs/yaw'
OUT = SOURCE / 'yaw_tracking_only'

def render_single_panel(summary, runs, cohort, figure_base=None):
    if figure_base is None:
        figure_base = SOURCE / 'yaw_absolute_error_ecdf_yaw_only'
    plt = renderer.plt
    fig, ax = plt.subplots(figsize=(6.2, 6.2*3122/3950))
    head_color = renderer.HEAD
    torso_color = renderer.TORSO
    fig.subplots_adjust(left=.13, right=.98, bottom=.13, top=.97)
    grid = summary.absoluteYawDeg.to_numpy()
    for segment, color, style in [('torso', torso_color, '-'),
                                   ('head', head_color, '-')]:
        curves = runs.pivot(index='absoluteYawDeg', columns='dataset',
                           values=segment+'PosteriorMedianTrajectoryCdf').reindex(index=grid, columns=cohort.dataset)
        ax.plot(grid, curves.to_numpy(), color=color, alpha=.22, lw=.7)
        mean = summary[segment+'RunEqualCdfMean'].to_numpy()
        se = summary[segment+'RunSE'].to_numpy()
        ax.fill_between(grid, np.maximum(0, mean-se), np.minimum(1, mean+se),
                        color=color, alpha=.13, lw=0)
        ax.plot(grid, mean, color=color, ls=style, lw=2.5)
    ax.set(xlim=(0, 180), ylim=(0, 1), xlabel='Absolute yaw error (deg)',
           ylabel='Cumulative proportion of frames', xticks=[0, 30, 60, 90, 120, 150, 180],
           yticks=np.linspace(0, 1, 11))
    ax.grid(alpha=.35, lw=.6)
    ax.tick_params(direction='in', labelsize=9)
    from matplotlib.ticker import StrMethodFormatter
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
    ax.axvline(10, color='.4', ls='--', lw=.8)
    ax.text(10.5, .005, '10 deg', rotation=90, va='bottom', ha='left', fontsize=9)
    handles = [renderer.Line2D([], [], color=head_color, lw=2.5, label='Head-lure yaw error'),
               renderer.Patch(facecolor=head_color, alpha=.13, label='Head: +/- 1 SE'),
               renderer.Line2D([], [], color=torso_color, lw=2.5, label='Torso-lure yaw error'),
               renderer.Patch(facecolor=torso_color, alpha=.13, label='Torso: +/- 1 SE')]
    ax.legend(handles=handles, loc='lower right', ncol=1, fontsize=8,
              frameon=True, fancybox=False, framealpha=1, facecolor='white',
              edgecolor='.25', handlelength=2.5, borderpad=.6)
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(figure_base.with_suffix('.'+extension), dpi=600, facecolor='white')
    plt.close(fig)
    (OUT / 'yaw_absolute_error_ecdf_caption.tex').write_text(
        'Absolute head--lure and torso--lure yaw errors across seven yaw-tracking pursuits '
        '(628 matched frames). Thick curves show equally weighted means of the per-run ECDFs; '
        'thin curves show individual runs and shading denotes one standard error across runs. '
        'Errors are absolute values of wrapped posterior-median yaw errors.\n')

def main():
    OUT.mkdir(exist_ok=True)
    cohort = pd.read_csv(SOURCE / 'yaw_absolute_error_cohort.csv')
    cohort = cohort[cohort.context.isin(['lateral pursuit', 'close pass-by pursuit'])]
    assert len(cohort) == 7
    frames = pd.read_csv(SOURCE / 'yaw_absolute_error_by_frame.csv')
    frames = frames[frames.dataset.isin(cohort.dataset)]
    assert len(frames) == 628
    assert frames.context.value_counts().to_dict() == {'lateral pursuit': 510, 'close pass-by pursuit': 118}
    runs = pd.read_csv(SOURCE / 'yaw_absolute_error_ecdf_by_run.csv')
    runs = runs[runs.dataset.isin(cohort.dataset)]
    grid = np.sort(runs.absoluteYawDeg.unique())
    summary = pd.DataFrame({'absoluteYawDeg': grid})
    for segment in ['head', 'torso']:
        curves = runs.pivot(index='absoluteYawDeg', columns='dataset',
                           values=segment+'PosteriorMedianTrajectoryCdf').reindex(index=grid, columns=cohort.dataset)
        assert curves.notna().all().all()
        for name in cohort.dataset:
            errors = frames.loc[frames.dataset == name, segment+'AbsoluteYawDeg'].to_numpy()
            independently_computed = (errors[:, None] <= grid[None, :]).mean(axis=0)
            assert np.allclose(curves[name], independently_computed, atol=1e-12)
        summary[segment+'RunEqualCdfMean'] = curves.mean(axis=1).to_numpy()
        summary[segment+'RunSE'] = (curves.std(axis=1, ddof=1)/np.sqrt(7)).to_numpy()
    for name, table in [('cohort', cohort), ('by_frame', frames), ('ecdf_by_run', runs), ('ecdf', summary)]:
        table.to_csv(OUT / f'yaw_absolute_error_{name}.csv', index=False)
    render_single_panel(summary, runs, cohort)
    (OUT / 'README.txt').write_text(
        'Seven yaw_tracking runs only: five lateral and two close pass-by pursuits; '
        f'628 matched frames, {cohort.animal.nunique()} cheetahs.\n'
        'Run-equal ECDF means and SE are recomputed for this subset. Each per-run curve '
        'was independently verified against the paired frame errors.\n'
        'Reproduce: python plotYawAbsoluteErrorEcdfYawOnly.py\n'
        'The original 18-run figure is preserved.\n')
    print(cohort[['dataset','context']].to_string(index=False))
    print(summary.loc[summary.absoluteYawDeg == 10].to_string(index=False))

if __name__ == '__main__':
    main()
