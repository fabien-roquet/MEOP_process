from __future__ import annotations

import numpy as np

import meop_process.processing.stabilise_sa_const_ct as stabilise_module


def test_stabilise_sp_const_ct_logs_and_skips_problem_profile(monkeypatch, capsys) -> None:
    max_float = np.finfo(float).max
    calls = {"count": 0}

    def fake_sigma0(sa, ct):
        _ = sa, ct
        calls["count"] += 1
        if calls["count"] == 1:
            return np.array([max_float, max_float * (1.0 - 1.0e-12)], dtype=float)
        return np.array([1.0, 2.0], dtype=float)

    monkeypatch.setattr(stabilise_module, "_sigma0_from_sa_ct", fake_sigma0)

    sp = np.array([[35.0, 35.0], [34.0, 34.0]], dtype=float)
    ct = np.array([[1.0, 1.0], [0.0, 0.0]], dtype=float)
    p = np.array([[10.0, 10.0], [20.0, 20.0]], dtype=float)

    out, metadata = stabilise_module.stabilise_SP_const_CT(
        sp,
        ct,
        p,
        return_metadata=True,
        smru_name="ct96-24-13",
    )

    np.testing.assert_allclose(out[:, 0], sp[:, 0])
    np.testing.assert_allclose(out[:, 1], sp[:, 1])
    assert metadata[0] is not None
    assert metadata[0].success is False
    assert metadata[0].profile_index == 0
    assert "skipped:" in metadata[0].status
    assert metadata[1] is not None
    assert metadata[1].success is True
    assert metadata[1].profile_index == 1

    captured = capsys.readouterr()
    assert "stabilisation skipped smru_name=ct96-24-13 profile_index=0:" in captured.out
