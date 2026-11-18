"""Check reproduced manuscript tables, preserving mismatches for inspection."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, default=ROOT / 'outputs')
    parser.add_argument('--include-spatial', action='store_true')
    args = parser.parse_args()
    frames = pd.read_csv(ROOT / 'reference/manuscript_frames.csv')
    pitch = pd.read_csv(args.output_root / 'pitch/tightness_run_equal/pitch_frame_level_errors.csv')
    pitch_keys = set(zip(pitch.datasetKey.str.split('/').str[-1], pitch.absFrame))
    pitch_ref = frames[frames.pitch_included]
    assert pitch_keys == set(zip(pitch_ref.dataset, pitch_ref.absolute_frame))
    assert len(pitch) == 1590 and pitch.datasetKey.nunique() == 16
    yaw = pd.read_csv(args.output_root / 'yaw/yaw_absolute_error_by_frame.csv')
    assert len(yaw) == 1708 and yaw.dataset.nunique() == 18
    assert set(zip(yaw.dataset,yaw.absFrame)) == set(zip(frames.dataset,frames.absolute_frame))
    assert yaw.groupby('context').size().to_dict() == {
        'frontal pursuit':1080, 'lateral pursuit':510, 'close pass-by pursuit':118}

    refs = ROOT / 'reference/analysis'
    for spec in json.loads((refs / 'checks.json').read_text()):
        if spec['optional'] and not args.include_spatial:
            continue
        actual = pd.read_csv(args.output_root / spec['output']).sort_values(spec['keys']).reset_index(drop=True)
        expected = pd.read_csv(refs / spec['reference']).sort_values(spec['keys']).reset_index(drop=True)
        assert len(actual) == len(expected), spec['output']
        pd.testing.assert_frame_equal(actual[spec['keys']], expected[spec['keys']], check_dtype=False)
        for column in expected.select_dtypes('number').columns:
            np.testing.assert_allclose(actual[column],expected[column],atol=1e-7,rtol=1e-7,
                equal_nan=True,err_msg=spec['output']+' :: '+column)
        print('PASS',spec['output'])

    states = pd.read_csv(args.output_root / 'yaw/geometry/gate_support_by_frame.csv')
    for gate, group in states.groupby('gateDeg'):
        assert gate in (5,10,15) and len(group) == 1708
        assert not group.duplicated(['dataset','absFrame']).any()
        p = group[['pAhead','pAligned','pLagged','pDirectionAmbiguous']].to_numpy()
        np.testing.assert_allclose(p.sum(axis=1),1,atol=1e-12)
        state = np.select((p[:,:3]>=.95).T,['ahead','aligned','lagged'],default='unresolved')
        assert np.array_equal(state,group.supportedState)
    assert set(states.gateDeg)=={5,10,15}
    print('PASS exact cohorts and posterior state classifications')


if __name__ == '__main__':
    main()
