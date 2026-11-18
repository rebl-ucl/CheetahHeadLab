#!/usr/bin/env python3
"""Equal-weight cheetahs after averaging their seven yaw-pursuit run ECDFs."""
from pathlib import Path
import numpy as np
import pandas as pd
import plotYawAbsoluteErrorEcdfYawOnly as plotter

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'outputs/yaw/yaw_tracking_only'
OUT = SOURCE / 'animal_equal'
STEM = 'yaw_absolute_error_ecdf_yaw_only_animal_equal'

def main():
    OUT.mkdir(exist_ok=True)
    cohort = pd.read_csv(SOURCE / 'yaw_absolute_error_cohort.csv')
    runs = pd.read_csv(SOURCE / 'yaw_absolute_error_ecdf_by_run.csv')
    assert len(cohort) == 7 and cohort.animal.nunique() == 4
    joined = runs.merge(cohort[['dataset','animal']], on='dataset', validate='many_to_one')
    columns = ['headPosteriorMedianTrajectoryCdf','torsoPosteriorMedianTrajectoryCdf']
    animals = joined.groupby(['animal','absoluteYawDeg'], as_index=False)[columns].mean()
    animals.to_csv(OUT / 'ecdf_by_animal.csv', index=False)
    grid = np.sort(animals.absoluteYawDeg.unique())
    summary = pd.DataFrame({'absoluteYawDeg':grid})
    for segment, column in zip(['head','torso'], columns):
        curves = animals.pivot(index='absoluteYawDeg',columns='animal',values=column).reindex(grid)
        assert curves.shape[1] == 4 and curves.notna().all().all()
        summary[segment+'AnimalEqualCdfMean'] = curves.mean(axis=1).to_numpy()
        summary[segment+'AnimalSE'] = curves.std(axis=1,ddof=1).to_numpy()/np.sqrt(4)
    summary.to_csv(OUT / 'ecdf_animal_equal.csv',index=False)
    # Reuse identical visual styling; pass animal curves as the plotting units.
    compatible = summary.rename(columns={s+a:s+b for s in ['head','torso'] for a,b in
        [('AnimalEqualCdfMean','RunEqualCdfMean'),('AnimalSE','RunSE')]})
    plotter.OUT = OUT
    plotter.render_single_panel(compatible, animals.rename(columns={'animal':'dataset'}),
                               pd.DataFrame({'dataset':cohort.animal.unique()}),
                               figure_base=ROOT.parent / 'outputs/yaw' / STEM)
    (OUT / 'yaw_absolute_error_ecdf_caption.tex').write_text(
        'Seven yaw pursuits, 628 frames, four cheetahs. Run ECDFs are averaged within '
        'each cheetah, then the four cheetahs receive equal weight. Thin curves show '
        'individual cheetahs; thick curves show means; shading denotes +/-1 SE across cheetahs.\n')
    print(summary.loc[summary.absoluteYawDeg==10].to_string(index=False))

if __name__ == '__main__':
    main()
