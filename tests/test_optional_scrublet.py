import pytest

import kimlabspatial.preprocessing as preprocessing


def test_missing_scrublet_only_blocks_doublet_detection(monkeypatch):
    monkeypatch.setattr(preprocessing, "scr", None)

    with pytest.raises(ImportError, match="Scrublet is required"):
        preprocessing.exclude_doublets(None, "", save=False)
