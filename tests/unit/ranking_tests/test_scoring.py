"""
Recommended testing commands:
    a) General:     pytest -v -p no:warnings --tb=no test_scoring.py
    b) Function:    pytest -v -p no:warnings --tb=no test_scoring.py::TestPScore::test_basic
    c) Doctest:     pytest -v --doctest-modules -p no:warnings test_scoring.py
"""
import pytest
import numpy as np

import xomics as xo


# Test p_score function
class TestPScore:
    def test_basic(self):
        # Test basic functionality
        result = xo.p_score(ids=['protein1', 'protein2'], x_fc=[2.4, 1.5], x_pvals=[0.05, 0.2])
        assert isinstance(result, np.ndarray)
        assert len(result) == 2

    def test_empty_input(self):
        # Test empty input
        with pytest.raises(ValueError):
            xo.p_score(ids=[], x_fc=[], x_pvals=[])

    def test_invalid_input_length(self):
        # Test different lengths of inputs
        with pytest.raises(ValueError):
            xo.p_score(ids=['protein1'], x_fc=[2.4, 1.5], x_pvals=[0.05])

    def test_non_numeric_input(self):
        # Test non-numeric input
        with pytest.raises(ValueError):
            xo.p_score(ids=['protein1', 'protein2'], x_fc=['a', 'b'], x_pvals=[0.05, 0.2])

    def test_negative_values(self):
        # Test negative fold change values
        result = xo.p_score(ids=['protein1', 'protein2'], x_fc=[-2.4, -1.5], x_pvals=[0.05, 0.2])
        assert np.allclose(result, [1.0, 0.], atol=1e-2)


# Test e_score function
class TestEScore:
    def test_basic(self):
        # Test basic functionality
        result = xo.e_score(
            ids=['protein1', 'protein2'],
            id_lists=[['protein1', 'protein2'], ['protein2']],
            x_fe=[2, 1.5],
            x_pvals=[0.05, 0.1]
        )
        assert isinstance(result, np.ndarray)
        assert len(result) == 2
        assert np.allclose(result, [0., 1.0], atol=1e-5)

    def test_empty_input(self):
        # Test empty input
        with pytest.raises(ValueError):
            xo.e_score(ids=[], id_lists=[], x_fe=[], x_pvals=[])

    def test_invalid_input_length(self):
        # Test different lengths of inputs
        with pytest.raises(ValueError):
            xo.e_score(ids=['protein1'], id_lists=[['protein1', 'protein2'], ['protein2']], x_fe=[2, 1.5], x_pvals=[0.05])

    def test_non_numeric_input(self):
        # Test non-numeric input
        with pytest.raises(ValueError):
            xo.e_score(
                ids=['protein1', 'protein2'],
                id_lists=[['protein1', 'protein2'], ['protein2']],
                x_fe=['a', 'b'],
                x_pvals=[0.05, 0.1]
            )

    def test_negative_values(self):
        # Test negative fold enrichment values
        with pytest.raises(ValueError):
            xo.e_score(
                ids=['protein1', 'protein2'],
                id_lists=[['protein1', 'protein2'], ['protein2']],
                x_fe=[-2, -1.5],
                x_pvals=[0.05, 0.1]
            )

    def test_invalid_ids(self):
        # Test for ids not in id_lists
        result = xo.e_score(
            ids=['protein3'],
            id_lists=[['protein1', 'protein2'], ['protein2']],
            x_fe=[2, 1.5],
            x_pvals=[0.05, 0.1]
        )
        assert result == [0]
