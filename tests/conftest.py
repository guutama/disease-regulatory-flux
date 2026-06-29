# Test-environment compatibility shim. ArviZ 0.17 (the newest for Python 3.9) imports
# scipy.signal.gaussian, which moved to scipy.signal.windows in SciPy >= 1.13. The pipeline
# pins a compatible SciPy at run time; this lets the test suite import ArviZ on a newer SciPy
# without changing it. It is a no-op when scipy.signal already provides gaussian.
try:
    import scipy.signal as _signal
    import scipy.signal.windows as _windows

    if not hasattr(_signal, "gaussian"):
        _signal.gaussian = _windows.gaussian
except Exception:
    pass
