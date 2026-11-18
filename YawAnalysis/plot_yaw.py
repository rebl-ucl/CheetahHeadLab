"""Render manuscript yaw panels after MATLAB run_yaw_analysis()."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotYawAbsoluteErrorEcdfYawOnly as ecdf
import plotYawAbsoluteErrorEcdfYawOnlyAnimalEqual as animal
import plotYawSignedErrorDistributionLateralOnly as signed
import plot_yaw_geometry as geometry

OUT = Path(__file__).resolve().parents[1] / 'outputs/yaw'
CONTEXTS = ['frontal pursuit', 'lateral pursuit', 'close pass-by pursuit']
STATES = ['ahead', 'aligned', 'lagged', 'unresolved']
COLORS = ['#0072B2', '#009E73', '#D55E00', '#AFAFAF']


def supported_modes():
    f = pd.read_csv(OUT / 'geometry/gate_support_by_frame.csv')
    expected = pd.read_csv(OUT / 'yaw_absolute_error_by_frame.csv')
    rows = []
    for gate, group in f.groupby('gateDeg'):
        assert len(group) == 1708
        assert set(zip(group.dataset, group.absFrame)) == set(zip(expected.dataset, expected.absFrame))
        p = group[['pAhead', 'pAligned', 'pLagged', 'pDirectionAmbiguous']].to_numpy()
        assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
        np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-12)
        classification = np.select((p[:, :3] >= .95).T, STATES[:3], default='unresolved')
        assert np.array_equal(classification, group.supportedState)
        for (name, context), run in group.groupby(['dataset', 'context'], sort=False):
            for state in STATES:
                rows.append(dict(dataset=name, context=context, gateDeg=gate,
                                 state=state, fraction=(run.supportedState == state).mean()))
    runs = pd.DataFrame(rows)
    summary = runs.groupby(['gateDeg', 'context', 'state'], sort=False).fraction.mean().mul(100).reset_index(name='percent')
    runs.to_csv(OUT / 'yaw_supported_modes_by_run.csv', index=False)
    summary.to_csv(OUT / 'yaw_supported_modes_summary.csv', index=False)

    def panel(ax, gate):
        v = summary[summary.gateDeg == gate].pivot(index='context', columns='state', values='percent').loc[CONTEXTS, STATES]
        np.testing.assert_allclose(v.sum(axis=1), 100)
        bottom = np.zeros(3)
        for state, color in zip(STATES, COLORS):
            h = v[state].to_numpy()
            ax.bar(range(3), h, bottom=bottom, width=.62, color=color, edgecolor='.2', linewidth=.7, label=state.capitalize())
            for i, value in enumerate(h):
                if value >= 3:
                    ax.text(i, bottom[i]+value/2, f'{value:.1f}%', ha='center', va='center', fontsize=10,
                            color='white' if state == 'ahead' else 'black')
            bottom += h
        ax.set(ylim=(0,100), ylabel='Run-equal proportion of frames (%)', xticks=range(3),
               xticklabels=['Frontal\npursuit', 'Lateral\npursuit', 'Close pass-by\npursuit'])
        ax.spines[['top','right']].set_visible(False)
    fig, ax = plt.subplots(figsize=(6.2,5.4), layout='constrained')
    panel(ax,10); ax.legend(loc='lower center', bbox_to_anchor=(.5,1), ncol=4, frameon=False, fontsize=9)
    for ext in ('png','pdf'):
        fig.savefig(OUT / f'yaw_mc_supported_modes_by_context.{ext}', dpi=400)
    plt.close(fig)
    fig, axes = plt.subplots(2,1,figsize=(6.2,10.8),layout='constrained')
    for ax, gate in zip(axes,[5,15]):
        panel(ax,gate); ax.set_title(f'{gate}° alignment gate', loc='left', fontsize=12)
    axes[0].legend(loc='lower center',bbox_to_anchor=(.5,1.12),ncol=4,frameon=False,fontsize=9)
    for ext in ('png','pdf'):
        fig.savefig(OUT / f'mode_threshold_sensitivity.{ext}', dpi=400)
    plt.close(fig)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    ecdf.main()
    animal.main()
    signed.main()
    geometry.main()
    supported_modes()
