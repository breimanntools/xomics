"""
This is testing script for the e_hit function.
"""
import pytest
import pandas as pd

# Import your actual p_score and e_score functions
from xomics import e_hits


# Test e_hits function
class TestEHits:
    def test_basic_functionality(self):
        result = e_hits(
            ids=['gene1', 'gene2', 'gene3'],
            id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']],
            terms=['term1', 'term2'],
            n_ids=2,
            n_terms=1
        )
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (1, 2)

    def test_empty_input(self):
        with pytest.raises(ValueError):
            e_hits(ids=[], id_lists=[], terms=[])

    #def test_invalid_id(self):
    #    with pytest.raises(ValueError):
    #        e_hits(ids=['gene4'], id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']], list_terms=['term1', 'term2'])

    def test_invalid_term_length(self):
        with pytest.raises(ValueError):
            e_hits(ids=['gene1', 'gene2', 'gene3'], id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']], terms=['term1'])

    def test_non_list_id_lists(self):
        with pytest.raises(ValueError):
            e_hits(ids=['gene1', 'gene2', 'gene3'], id_lists=['gene1', 'gene2', 'gene3'], terms=['term1', 'term2'])

    def test_valid_n_ids(self):
        result = e_hits(
            ids=['gene1', 'gene2', 'gene3'],
            id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']],
            terms=['term1', 'term2'],
            n_ids=2,
            n_terms=None
        )
        assert result.shape == (2, 2)

    def test_valid_n_terms(self):
        result = e_hits(
            ids=['gene1', 'gene2', 'gene3'],
            id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']],
            terms=['term1', 'term2'],
            n_ids=None,
            n_terms=1
        )
        assert result.shape == (1, 3)

    def test_valid_n_ids_and_n_terms(self):
        result = e_hits(
            ids=['gene1', 'gene2', 'gene3'],
            id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']],
            terms=['term1', 'term2'],
            n_ids=2,
            n_terms=2
        )
        assert result.shape == (2, 2)

    def test_invalid_n_ids(self):
        with pytest.raises(ValueError):
            e_hits(ids=['gene1', 'gene2', 'gene3'], id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']], terms=['term1', 'term2'], n_ids=-1)

    def test_invalid_n_terms(self):
        with pytest.raises(ValueError):
            e_hits(ids=['gene1', 'gene2', 'gene3'], id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']], terms=['term1', 'term2'], n_terms=-1)
