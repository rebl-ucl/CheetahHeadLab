"""Check that the full analysis workflow exports only inventoried manuscript panels."""
from pathlib import Path
import argparse
import json

ROOT = Path(__file__).resolve().parents[1]


def check(output_root):
    inventory = json.loads((ROOT / 'reference/manuscript_figures.json').read_text())
    panels = inventory['panels']
    expected = {p['stem'] for p in panels}
    assert len(expected) == len(panels), 'Duplicate panel paths in figure inventory'
    missing = sorted(stem + '.png' for stem in expected
                     if not (output_root / (stem + '.png')).is_file()
                     or (output_root / (stem + '.png')).stat().st_size == 0)
    extras = []
    for folder in ('pitch', 'yaw'):
        for path in (output_root / folder).rglob('*'):
            if path.is_file() and path.suffix.lower() in inventory['allowed_formats']:
                stem = path.relative_to(output_root).with_suffix('').as_posix()
                if stem not in expected:
                    extras.append(path.relative_to(output_root).as_posix())
    if missing or extras:
        parts = []
        if missing:
            parts.append('Missing manuscript panels:\n  ' + '\n  '.join(missing))
        if extras:
            parts.append('Non-manuscript exports (including possible leftovers):\n  '
                         + '\n  '.join(sorted(extras)))
        raise SystemExit('\n'.join(parts))
    print(f'PASS: {len(expected)} manuscript panels, no unlisted figure exports.')
    print('Composite illustrations and uncovered panels: docs/FIGURE_COVERAGE.md')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, default=ROOT / 'outputs')
    check(parser.parse_args().output_root)
