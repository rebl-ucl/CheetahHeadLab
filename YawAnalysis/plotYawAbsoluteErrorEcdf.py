#!/usr/bin/env python3
"""Render the current verified ECDF tables without headless MATLAB font issues."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parents[1]/'outputs'/'yaw'
HEAD='#0072B2'
TORSO='#E69F00'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                    'axes.spines.top':False,'axes.spines.right':False,
                    'svg.fonttype':'none','pdf.fonttype':42})

