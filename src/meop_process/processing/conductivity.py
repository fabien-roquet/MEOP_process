from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

try:  # pragma: no cover - exercised only when gsw is installed
    import gsw as _gsw  # type: ignore
except Exception:  # pragma: no cover - current container does not provide gsw
    _gsw = None


GSW_SSO = 35.16504
PSS78_C35 = 42.9140
THERMAL_ALPHA = (0.09, 0.05)
THERMAL_BETA = 0.06
THERMAL_DT = 1.0
SMOOTHING_KERNEL = np.array([1.0, 4.0, 6.0, 4.0, 1.0], dtype=np.float64) / 16.0


@dataclass(frozen=True)
class ThermalCorrectionResult:
    temperature: np.ndarray
    salinity: np.ndarray
    highpass_temperature: np.ndarray
    method: str


def _broadcast_to_shape(values: np.ndarray | float, shape: tuple[int, ...]) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.shape == shape:
        return array
    return np.broadcast_to(array, shape).astype(np.float64)


def sp_from_c(c: np.ndarray, t: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Practical salinity from conductivity using the PSS-78 expression.

    This follows the core branch of ``gsw_SP_from_C``. The low-salinity Hill correction is
    intentionally omitted because MEOP processing thresholds already operate in the open-ocean
    salinity range where SP < 2 is not expected.
    """

    conductivity = np.asarray(c, dtype=np.float64)
    temp = _broadcast_to_shape(np.asarray(t, dtype=np.float64), conductivity.shape)
    pres = _broadcast_to_shape(np.asarray(p, dtype=np.float64), conductivity.shape)

    out = np.full(conductivity.shape, np.nan, dtype=np.float64)
    mask = np.isfinite(conductivity) & np.isfinite(temp) & np.isfinite(pres) & (conductivity > 0)
    if not np.any(mask):
        return out

    a0, a1, a2, a3, a4, a5 = 0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081
    b0, b1, b2, b3, b4, b5 = 0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144
    c0, c1, c2, c3, c4 = 0.6766097, 2.00564e-2, 1.104259e-4, -6.9698e-7, 1.0031e-9
    d1, d2, d3, d4 = 3.426e-2, 4.464e-4, 4.215e-1, -3.107e-3
    e1, e2, e3 = 2.070e-5, -6.370e-10, 3.989e-15
    k = 0.0162

    c_use = conductivity[mask]
    t_use = temp[mask]
    p_use = pres[mask]

    t68 = t_use * 1.00024
    ft68 = (t68 - 15.0) / (1.0 + k * (t68 - 15.0))
    r = c_use / PSS78_C35
    rt_lc = c0 + (c1 + (c2 + (c3 + c4 * t68) * t68) * t68) * t68
    rp = 1.0 + (p_use * (e1 + e2 * p_use + e3 * p_use * p_use)) / (
        1.0 + d1 * t68 + d2 * t68 * t68 + (d3 + d4 * t68) * r
    )
    rt = r / (rp * rt_lc)
    rt = np.where(rt > 0.0, rt, np.nan)
    rtx = np.sqrt(rt)

    sp = a0 + (a1 + (a2 + (a3 + (a4 + a5 * rtx) * rtx) * rtx) * rtx) * rtx
    sp = sp + ft68 * (b0 + (b1 + (b2 + (b3 + (b4 + b5 * rtx) * rtx) * rtx) * rtx) * rtx)
    sp = np.where(sp < 0.0, 0.0, sp)
    out[mask] = sp
    return out


def c_from_sp(sp: np.ndarray, t: np.ndarray, p: np.ndarray, *, iterations: int = 6) -> np.ndarray:
    """Conductivity from practical salinity.

    Prefer the Python ``gsw`` implementation when available. Otherwise solve the PSS-78
    inversion numerically against :func:`sp_from_c`.
    """

    practical_salinity = np.asarray(sp, dtype=np.float64)
    temp = _broadcast_to_shape(np.asarray(t, dtype=np.float64), practical_salinity.shape)
    pres = _broadcast_to_shape(np.asarray(p, dtype=np.float64), practical_salinity.shape)

    if _gsw is not None:  # pragma: no branch - deterministic branch when dependency is installed
        try:
            return np.asarray(_gsw.C_from_SP(practical_salinity, temp, pres), dtype=np.float64)
        except Exception:
            pass

    out = np.full(practical_salinity.shape, np.nan, dtype=np.float64)
    mask = np.isfinite(practical_salinity) & np.isfinite(temp) & np.isfinite(pres) & (practical_salinity >= 0.0)
    if not np.any(mask):
        return out

    target = practical_salinity[mask]
    t_use = temp[mask]
    p_use = pres[mask]
    guess = PSS78_C35 * np.clip(target, 0.0, 42.0) / 35.0
    guess = np.where(target == 0.0, 0.0, np.maximum(guess, 1e-6))

    for _ in range(iterations):
        estimate = sp_from_c(guess, t_use, p_use)
        delta = np.where(np.isfinite(estimate), estimate - target, 0.0)
        step = np.maximum(1e-3, np.abs(guess) * 1e-5)
        d_plus = sp_from_c(guess + step, t_use, p_use)
        d_minus = sp_from_c(np.maximum(guess - step, 1e-9), t_use, p_use)
        deriv = (d_plus - d_minus) / (2.0 * step)
        valid = np.isfinite(delta) & np.isfinite(deriv) & (np.abs(deriv) > 1e-10)
        if not np.any(valid):
            break
        guess[valid] = guess[valid] - delta[valid] / deriv[valid]
        guess = np.where(np.isfinite(guess), np.maximum(guess, 0.0), guess)

    out[mask] = guess
    return out


def smooth_profile(values: np.ndarray) -> np.ndarray:
    padded = np.concatenate([np.full(2, np.nan), values.astype(np.float64), np.full(2, np.nan)])
    smoothed = np.full(values.shape, np.nan, dtype=np.float64)
    for idx in range(values.size):
        window = padded[idx : idx + 5]
        if np.all(np.isfinite(window)):
            smoothed[idx] = float(np.dot(SMOOTHING_KERNEL, window))
    replace = np.isnan(smoothed) & np.isfinite(values)
    smoothed[replace] = values[replace]
    return smoothed


def _thermal_cell_profile(
    sp: np.ndarray,
    temp: np.ndarray,
    pres: np.ndarray,
    *,
    alpha: tuple[float, float] = THERMAL_ALPHA,
    beta: float = THERMAL_BETA,
    dt: float = THERMAL_DT,
    direction: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if np.count_nonzero(np.isfinite(sp) & np.isfinite(temp) & np.isfinite(pres)) < 3:
        return temp.copy(), sp.copy(), np.full(temp.shape, np.nan, dtype=np.float64)

    work_p = pres[::-1].copy() if direction else pres.copy()
    work_t = temp[::-1].copy() if direction else temp.copy()
    work_s = sp[::-1].copy() if direction else sp.copy()

    conductivity = c_from_sp(work_s, work_t, work_p)
    tau = 1.0 / beta
    response = tau / (tau + dt)

    thp = np.full(work_t.shape, np.nan, dtype=np.float64)
    thp[0] = 0.0 if np.isfinite(work_t[0]) else np.nan
    for idx in range(1, work_t.size):
        if not np.isfinite(work_t[idx - 1]):
            thp[idx] = 0.0
            continue
        if not np.isfinite(work_t[idx]):
            thp[idx] = np.nan
            continue
        previous = 0.0 if not np.isfinite(thp[idx - 1]) else thp[idx - 1]
        thp[idx] = response * (previous + work_t[idx] - work_t[idx - 1])

    grad_c = c_from_sp(work_s, work_t + 0.5, work_p) - c_from_sp(work_s, work_t - 0.5, work_p)
    grad_c = np.where(np.isfinite(grad_c), grad_c, 0.0)

    factor = 1.0 / (1.0 - 0.5 * beta * dt)
    tc = alpha[0] * factor * thp
    cc = alpha[1] * factor * grad_c * thp

    t_corr = work_t + tc
    c_corr = conductivity + cc
    s_corr = sp_from_c(c_corr, t_corr, work_p)

    if direction:
        return t_corr[::-1], s_corr[::-1], thp[::-1]
    return t_corr, s_corr, thp


def thermal_cell_correction(
    sp: np.ndarray,
    temp: np.ndarray,
    pres: np.ndarray,
    *,
    alpha: tuple[float, float] = THERMAL_ALPHA,
    beta: float = THERMAL_BETA,
    dt: float = THERMAL_DT,
    direction: bool = True,
) -> ThermalCorrectionResult:
    practical_salinity = np.asarray(sp, dtype=np.float64)
    in_situ_temp = np.asarray(temp, dtype=np.float64)
    pressure = np.asarray(pres, dtype=np.float64)
    if practical_salinity.shape != in_situ_temp.shape or practical_salinity.shape != pressure.shape:
        raise ValueError("sp, temp, and pres must have the same shape")

    t_corr = np.full(practical_salinity.shape, np.nan, dtype=np.float64)
    s_corr = np.full(practical_salinity.shape, np.nan, dtype=np.float64)
    thp = np.full(practical_salinity.shape, np.nan, dtype=np.float64)

    for profile_index in range(practical_salinity.shape[0]):
        corrected_t, corrected_s, hp = _thermal_cell_profile(
            practical_salinity[profile_index],
            in_situ_temp[profile_index],
            pressure[profile_index],
            alpha=alpha,
            beta=beta,
            dt=dt,
            direction=direction,
        )
        t_corr[profile_index, :] = corrected_t
        s_corr[profile_index, :] = corrected_s
        thp[profile_index, :] = hp

    method = "gsw" if _gsw is not None else "pss78-fallback"
    return ThermalCorrectionResult(temperature=t_corr, salinity=s_corr, highpass_temperature=thp, method=method)
