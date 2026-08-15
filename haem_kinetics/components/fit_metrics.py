"""Fit diagnostics vs Garnie-style heme fractionation (fg Fe/cell).

These scores track whether a mechanistic change moves the *stated* ODEs
toward the assay. They are not an objective to minimize by adding fudges.
"""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional

import numpy as np
import pandas as pd

from haem_kinetics.components.experimental_data import ExperimentalData


DEFAULT_SPECIES_MAP = {
    'Hb': 'conc_hb_dv',
    'Hm': 'conc_fe3pp',
    'Hz': 'conc_hz',
}

HOST_KEY = 'conc_hb_rbc'
SKIP_TOTAL = {'conc_fe3pp_free', 'conc_hb_assay', 'conc_hb_dv_obs'}


def _interp_at(series: pd.Series, times: np.ndarray) -> np.ndarray:
    s = series.sort_index()
    x = s.index.to_numpy(dtype=float)
    y = s.to_numpy(dtype=float)
    out = np.full(len(times), np.nan, dtype=float)
    inside = (times >= x[0]) & (times <= x[-1])
    out[inside] = np.interp(times[inside], x, y)
    return out


def _series_metrics(pred: np.ndarray, obs: np.ndarray, sem: np.ndarray) -> Dict[str, float]:
    mask = np.isfinite(pred) & np.isfinite(obs)
    n = int(mask.sum())
    if n == 0:
        return {'n': 0, 'rmse': float('nan'), 'mae': float('nan'),
                'mean_signed_error': float('nan'), 'chi2_red': float('nan')}
    err = pred[mask] - obs[mask]
    rmse = float(np.sqrt(np.mean(err ** 2)))
    mae = float(np.mean(np.abs(err)))
    mean_signed_error = float(np.mean(err))
    sem_ok = mask & np.isfinite(sem) & (sem > 0.0)
    n_chi = int(sem_ok.sum())
    if n_chi == 0:
        chi2_red = float('nan')
    else:
        chi2_red = float(np.mean(((pred[sem_ok] - obs[sem_ok]) / sem[sem_ok]) ** 2))
    return {'n': n, 'rmse': rmse, 'mae': mae,
            'mean_signed_error': mean_signed_error, 'chi2_red': chi2_red}


def dv_total_fe_fg(df: pd.DataFrame) -> pd.Series:
    cols = [c for c in df.columns
            if c.startswith('conc_') and c not in SKIP_TOTAL and c != HOST_KEY]
    return df[cols].sum(axis=1)


def score_fractionation(
    model_fg: pd.DataFrame,
    exp_data: ExperimentalData,
    species_map: Optional[Mapping[str, str]] = None,
    free_haem_cols: Optional[List[str]] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Interpolate the model onto experimental ages and score Hb, Hm, Hz, DV Fe.

    RMSE/MAE/mean signed error are fg/cell. chi2_red is mean(((pred−obs)/SEM)²);
    ~1 means residuals match reported assay scatter. mean signed error
    = mean(pred − obs) (not mean squared error); > 0 means the model is high.
    """
    if exp_data is None or exp_data.data.empty:
        raise ValueError('exp_data has no fractionation table')
    if model_fg.empty:
        raise ValueError('model concentrations are empty; run() first')

    mapping = dict(species_map or DEFAULT_SPECIES_MAP)
    df = model_fg.copy()
    if free_haem_cols:
        df['conc_fe3pp_free'] = df[free_haem_cols].sum(axis=1)
        mapping['Hm'] = 'conc_fe3pp_free'

    times = exp_data.data.index.to_numpy(dtype=float)
    out: Dict[str, Dict[str, float]] = {}

    for exp_col, model_col in mapping.items():
        if model_col not in df.columns:
            continue
        obs = exp_data.data[exp_col].to_numpy(dtype=float)
        sem = exp_data.data[f'{exp_col}:SEM'].to_numpy(dtype=float)
        pred = _interp_at(df[model_col], times)
        out[exp_col] = _series_metrics(pred, obs, sem)

    exp_dv = (
        exp_data.data['Hb'].to_numpy(dtype=float)
        + exp_data.data['Hm'].to_numpy(dtype=float)
        + exp_data.data['Hz'].to_numpy(dtype=float)
    )
    exp_dv_sem = np.sqrt(
        exp_data.data['Hb:SEM'].to_numpy(dtype=float) ** 2
        + exp_data.data['Hm:SEM'].to_numpy(dtype=float) ** 2
        + exp_data.data['Hz:SEM'].to_numpy(dtype=float) ** 2
    )
    pred_dv = _interp_at(dv_total_fe_fg(df), times)
    out['DV_Fe'] = _series_metrics(pred_dv, exp_dv, exp_dv_sem)
    return out


def format_fit_metrics(metrics: Mapping[str, Mapping[str, float]], indent: str = '  ') -> str:
    lines = [f'{indent}species   RMSE   MAE  signed  chi2_red   n']
    for name in ('Hb', 'Hm', 'Hz', 'DV_Fe'):
        if name not in metrics:
            continue
        m = metrics[name]
        lines.append(
            f'{indent}{name:<8} {m["rmse"]:6.2f} {m["mae"]:5.2f} '
            f'{m["mean_signed_error"]:7.2f} {m["chi2_red"]:8.2f} {int(m["n"]):3d}'
        )
    return '\n'.join(lines)
