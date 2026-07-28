"""Installation smoke test for the public notebook environment."""


def test_core_notebook_imports():
    import anndata  # noqa: F401
    import dask  # noqa: F401
    import numpy  # noqa: F401
    import pandas  # noqa: F401
    import scanpy  # noqa: F401
    import scipy  # noqa: F401
    import spatialdata  # noqa: F401
    import squidpy  # noqa: F401
    import zarr  # noqa: F401

    import kimlabspatial.preprocessing  # noqa: F401
    import scale_aware_st  # noqa: F401
