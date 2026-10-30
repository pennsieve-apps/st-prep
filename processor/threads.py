"""
Thread count from THREADS (0 = the cgroup CPU quota, never the host's core count), set as the
BLAS / OpenMP / numba environment variables. Import this before numpy.
"""
from __future__ import annotations

import math
import os

VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS")


def cgroup_cpus() -> int | None:
    """CPUs from the cgroup quota (v2 cpu.max, v1 cfs_quota_us / cfs_period_us), or None without a quota."""
    try:
        with open("/sys/fs/cgroup/cpu.max") as fh:
            quota, period = fh.read().split()[:2]
        if quota != "max":
            return max(1, math.ceil(int(quota) / int(period)))
        return None
    except (OSError, ValueError):
        pass
    try:
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as fh:
            quota = int(fh.read())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as fh:
            period = int(fh.read())
        if quota > 0 and period > 0:
            return max(1, math.ceil(quota / period))
    except (OSError, ValueError):
        pass
    return None


def configure(value: str | None = None) -> int:
    """Set the thread environment variables from THREADS (or the quota; 2 when there is none). Returns the count."""
    raw = (value if value is not None else os.environ.get("THREADS", "0")) or "0"
    try:
        n = int(raw)
    except ValueError:
        n = 0
    if n <= 0:
        n = cgroup_cpus() or 2
    for v in VARS:
        os.environ[v] = str(n)
    return n


THREADS = configure()
