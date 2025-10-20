import os
from threadpoolctl import threadpool_limits

# Environment-level fallback (before any import)
for var in ["OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
            "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(var, "1")

# Runtime enforcement (in case BLAS already initialized)
threadpool_limits(limits=1, user_api='blas')

