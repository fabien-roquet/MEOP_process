"""Stabilise a water column by minimally adjusting salinity while keeping CT fixed.

This module provides a direct Python port of the optimization-based TEOS-10 workflow when the
``gsw`` package is available, plus a lightweight fallback that keeps the development and test
workflow usable in offline environments.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

try:  # pragma: no cover - exercised only when gsw is installed
    import gsw as _gsw  # type: ignore
except Exception:  # pragma: no cover - current container does not provide gsw
    _gsw = None

from .conductivity import GSW_SSO
from .qc import _sigma0_eos80_fallback

ArrayLike = np.ndarray | list[float] | tuple[float, ...]
SA_SCALE = GSW_SSO / 35.0


@dataclass(frozen=True)
class SolverInfo:
    solver: str
    status: str
    success: bool


class StabilisationError(RuntimeError):
    """Raised when stabilisation fails."""


def _as_2d_column_major(x: ArrayLike, name: str) -> tuple[np.ndarray, bool]:
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 1:
        return arr[:, None], True
    if arr.ndim == 2:
        return arr, False
    raise ValueError(f"{name} must be 1D or 2D, got ndim={arr.ndim}")


def _broadcast_p(p: ArrayLike, shape: tuple[int, int]) -> np.ndarray:
    m, n = shape
    p_arr = np.asarray(p, dtype=float)
    if p_arr.ndim == 0:
        raise ValueError("p must not be scalar")
    if p_arr.ndim == 1:
        if p_arr.size == m:
            return np.repeat(p_arr[:, None], n, axis=1)
        if p_arr.size == n:
            return np.repeat(p_arr[None, :], m, axis=0)
        raise ValueError(f"1D p must have length M={m} or N={n}, got length {p_arr.size}")
    if p_arr.ndim == 2:
        if p_arr.shape == (m, n):
            return p_arr
        if p_arr.shape == (m, 1):
            return np.repeat(p_arr, n, axis=1)
        if p_arr.shape == (1, n):
            return np.repeat(p_arr, m, axis=0)
        if p_arr.shape == (n, 1):
            return np.repeat(p_arr.T, m, axis=0)
        if p_arr.shape == (1, m):
            return np.repeat(p_arr.T, n, axis=1)
    raise ValueError(f"p has incompatible shape {p_arr.shape}; expected broadcastable to {(m, n)}")


def _broadcast_n2_lowerlimit(n2_lowerlimit: float | ArrayLike, shape: tuple[int, int]) -> np.ndarray:
    m, n = shape
    k = m - 1
    if np.isscalar(n2_lowerlimit):
        return np.full((k, n), float(n2_lowerlimit), dtype=float)
    arr = np.asarray(n2_lowerlimit, dtype=float)
    if arr.ndim == 1:
        if arr.size == k:
            return np.repeat(arr[:, None], n, axis=1)
        if arr.size == n:
            return np.repeat(arr[None, :], k, axis=0)
        raise ValueError(f"1D n2_lowerlimit must have length M-1={k} or N={n}, got {arr.size}")
    if arr.ndim == 2:
        if arr.shape == (k, n):
            return arr
        if arr.shape == (k, 1):
            return np.repeat(arr, n, axis=1)
        if arr.shape == (1, n):
            return np.repeat(arr, k, axis=0)
    raise ValueError(f"n2_lowerlimit has incompatible shape {arr.shape}")


def _first_difference_matrix(n: int) -> np.ndarray:
    a = np.zeros((n - 1, n), dtype=float)
    rows = np.arange(n - 1)
    a[rows, rows] = 1.0
    a[rows, rows + 1] = -1.0
    return a


def _select_valid_profile(sa: np.ndarray, ct: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mask = np.isfinite(sa) & np.isfinite(ct) & np.isfinite(p)
    idx = np.flatnonzero(mask)
    return idx, sa[idx], ct[idx], p[idx]


def nsquared_min_const_CT(
    sa: ArrayLike,
    ct: ArrayLike,
    p: ArrayLike,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sa_arr = np.asarray(sa, dtype=float)
    ct_arr = np.asarray(ct, dtype=float)
    p_arr = np.asarray(p, dtype=float)
    if sa_arr.ndim != 1 or ct_arr.ndim != 1 or p_arr.ndim != 1:
        raise ValueError("SA, CT, and p must be 1D arrays")
    if not (sa_arr.size == ct_arr.size == p_arr.size):
        raise ValueError("SA, CT, and p must have the same length")
    if sa_arr.size < 2:
        raise ValueError("At least 2 levels are required")
    dp = np.diff(p_arr)
    if np.any(dp <= 0):
        raise ValueError("Pressure must be strictly increasing")
    dsa = np.diff(sa_arr)
    dct = np.diff(ct_arr)
    if _gsw is None:
        raise StabilisationError("gsw is required for the TEOS-10 stabilisation solver")
    specvol, alpha, beta = _gsw.specvol_alpha_beta(sa_arr, ct_arr, p_arr)
    _ = specvol
    grav2 = 9.7963 ** 2
    db2pa = 1e4
    rho0 = 1e3
    prefactor = rho0 * grav2 / (db2pa * dp)
    n2_shallow = prefactor * (beta[:-1] * dsa - alpha[:-1] * dct)
    n2_deep = prefactor * (beta[1:] * dsa - alpha[1:] * dct)
    use_deep = n2_deep < n2_shallow
    n2 = np.where(use_deep, n2_deep, n2_shallow)
    n2_alpha = np.where(use_deep, alpha[1:], alpha[:-1])
    n2_beta = np.where(use_deep, beta[1:], beta[:-1])
    return n2, n2_alpha, n2_beta, dsa, dct, dp


def _solve_qp_scipy(a: np.ndarray, b_u: np.ndarray, x_lower: np.ndarray) -> tuple[np.ndarray, SolverInfo]:
    from scipy.optimize import Bounds, LinearConstraint, minimize

    n = x_lower.size

    def fun(x: np.ndarray) -> float:
        return 0.5 * float(np.dot(x, x))

    def jac(x: np.ndarray) -> np.ndarray:
        return x

    def hess(x: np.ndarray) -> np.ndarray:
        return np.eye(n)

    linear_constraint = LinearConstraint(a, -np.inf * np.ones_like(b_u), b_u)
    bounds = Bounds(lb=x_lower, ub=np.full(n, np.inf))
    x0 = np.maximum(np.zeros(n, dtype=float), x_lower)
    res = minimize(
        fun=fun,
        x0=x0,
        method="trust-constr",
        jac=jac,
        hess=hess,
        constraints=[linear_constraint],
        bounds=bounds,
        options={"verbose": 0},
    )
    if not res.success or res.x is None:
        raise StabilisationError(f"SciPy trust-constr failed: {res.message}")
    return np.asarray(res.x, dtype=float), SolverInfo(solver="scipy-trust-constr", status=str(res.message), success=True)


def _pav_non_decreasing(values: np.ndarray) -> np.ndarray:
    y = values.astype(float)
    n = y.size
    solution = y.copy()
    weights = np.ones(n, dtype=float)
    i = 0
    while i < n - 1:
        if solution[i] <= solution[i + 1]:
            i += 1
            continue
        total = solution[i] * weights[i] + solution[i + 1] * weights[i + 1]
        weight = weights[i] + weights[i + 1]
        solution[i] = solution[i + 1] = total / weight
        weights[i] = weights[i + 1] = weight
        j = i
        while j > 0 and solution[j - 1] > solution[j]:
            total = solution[j - 1] * weights[j - 1] + solution[j] * weights[j]
            weight = weights[j - 1] + weights[j]
            solution[j - 1] = solution[j] = total / weight
            weights[j - 1] = weights[j] = weight
            j -= 1
        i += 1
    return np.maximum.accumulate(solution)


def _sigma0_from_sa_ct(sa: np.ndarray, ct: np.ndarray) -> np.ndarray:
    if _gsw is not None:
        return np.asarray(_gsw.sigma0(sa, ct), dtype=float)
    sp = np.asarray(sa, dtype=float) / SA_SCALE
    return _sigma0_eos80_fallback(sp, ct)


def _invert_sigma0_to_salinity(target_sigma: np.ndarray, ct: np.ndarray, initial_sa: np.ndarray) -> np.ndarray:
    sa = initial_sa.astype(float).copy()
    for _ in range(10):
        sigma = _sigma0_from_sa_ct(sa, ct)
        eps = 1e-3
        sigma_plus = _sigma0_from_sa_ct(sa + eps, ct)
        sigma_minus = _sigma0_from_sa_ct(np.maximum(sa - eps, 0.0), ct)
        deriv = (sigma_plus - sigma_minus) / (2.0 * eps)
        valid = np.isfinite(target_sigma) & np.isfinite(ct) & np.isfinite(sa) & np.isfinite(deriv) & (np.abs(deriv) > 1e-10)
        if not np.any(valid):
            break
        sa[valid] = sa[valid] - (sigma[valid] - target_sigma[valid]) / deriv[valid]
        sa = np.where(np.isfinite(sa), np.maximum(sa, 0.0), sa)
    return sa


def _stabilise_sa_isotonic(sa_in: np.ndarray, ct: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, SolverInfo]:
    _ = p
    sigma = _sigma0_from_sa_ct(sa_in, ct)
    if not np.isfinite(sigma).any() or np.all(np.diff(sigma[np.isfinite(sigma)]) >= 0.0):
        solver = "isotonic-gsw" if _gsw is not None else "isotonic-eos80"
        return sa_in.copy(), SolverInfo(solver=solver, status="already_stable", success=True)
    target = _pav_non_decreasing(sigma)
    target += np.linspace(0.0, 1e-7, target.size)
    sa_stable = _invert_sigma0_to_salinity(target, ct, sa_in)
    solver = "isotonic-gsw" if _gsw is not None else "isotonic-eos80"
    return sa_stable, SolverInfo(solver=solver, status="approximate", success=True)


def stabilise_SA_const_CT(
    SA_in: ArrayLike,
    CT: ArrayLike,
    p: ArrayLike,
    n2_lowerlimit: float | ArrayLike = 1e-9,
    solver: Literal["auto", "scipy", "isotonic"] = "auto",
    copy: bool = True,
    return_metadata: bool = False,
) -> np.ndarray | tuple[np.ndarray, list[SolverInfo | None]]:
    sa, sa_was_1d = _as_2d_column_major(SA_in, "SA_in")
    ct_arr, ct_was_1d = _as_2d_column_major(CT, "CT")
    if sa.shape != ct_arr.shape:
        raise ValueError(f"SA_in and CT must have the same shape, got {sa.shape} and {ct_arr.shape}")
    m, n = sa.shape
    if m < 2:
        raise ValueError("At least 2 vertical levels are required")
    p_arr = _broadcast_p(p, (m, n))
    n2_min_arr = _broadcast_n2_lowerlimit(n2_lowerlimit, (m, n))
    sa_out = sa.copy() if copy else sa
    metadata: list[SolverInfo | None] = [None] * n
    c = 1.250423402612047e02
    rho0 = 1e3

    for j in range(n):
        idx, sa_prof, ct_prof, p_prof = _select_valid_profile(sa[:, j], ct_arr[:, j], p_arr[:, j])
        if sa_prof.size < 2:
            continue
        if solver in ("auto", "isotonic"):
            sa_stable, info = _stabilise_sa_isotonic(sa_prof, ct_prof, p_prof)
            sa_out[idx, j] = sa_stable
            metadata[j] = info
            continue
        if _gsw is None:
            raise StabilisationError("gsw is required for scipy TEOS-10 stabilisation")
        n2, n2_alpha, n2_beta, dsa, dct, dp = nsquared_min_const_CT(sa_prof, ct_prof, p_prof)
        n2_min_prof_full = n2_min_arr[:, j]
        n2_min_prof = n2_min_prof_full[idx[1:] - 1]
        unstable = (n2 - n2_min_prof) < 0.0
        if not np.any(unstable):
            continue
        b_u = dsa - (n2_alpha / n2_beta) * dct - (c / rho0) * (n2_min_prof * dp / n2_beta)
        pl = sa_prof.size
        a = _first_difference_matrix(pl)
        x_lower = -sa_prof
        x, info = _solve_qp_scipy(a=a, b_u=b_u, x_lower=x_lower)
        sa_out[idx, j] = sa_prof + x
        metadata[j] = info

    if sa_was_1d and ct_was_1d:
        sa_out = sa_out[:, 0]
    if return_metadata:
        return sa_out, metadata
    return sa_out


def stabilise_SP_const_CT(
    SP_in: ArrayLike,
    CT: ArrayLike,
    p: ArrayLike,
    **kwargs,
) -> np.ndarray | tuple[np.ndarray, list[SolverInfo | None]]:
    sp = np.asarray(SP_in, dtype=float)
    sa = sp * SA_SCALE
    result = stabilise_SA_const_CT(sa, CT, p, **kwargs)
    if isinstance(result, tuple):
        sa_out, metadata = result
        return np.asarray(sa_out, dtype=float) / SA_SCALE, metadata
    return np.asarray(result, dtype=float) / SA_SCALE
