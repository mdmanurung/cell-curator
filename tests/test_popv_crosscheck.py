from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import cell_curator.crosschecks as crosschecks
from cell_curator.contracts import ContractError
from cell_curator.crosschecks import celltypist_hypothesis, popv_hypothesis, predict_popv_file
from cell_curator.reconciliation import reconcile_predictions


def _input() -> ad.AnnData:
    value = ad.AnnData(
        np.array([[30, 3], [2, 31], [1, 1]], dtype=np.float32),
        obs=pd.DataFrame({'cluster': ['a', 'b', 'a']}, index=['c0', 'c1', 'c2']),
        var=pd.DataFrame(index=['G1', 'G2']),
    )
    value.layers['counts'] = value.X.copy()
    return value


def _bundle(path: Path) -> Path:
    path.mkdir()
    (path / 'preprocessing.json').write_text(json.dumps({'gene_names': ['G1', 'G2']}))
    (path / 'scvi').mkdir()
    (path / 'scvi/model.pt').write_bytes(b'local scvi bundle')
    for name in ['celltypist.pkl', 'svm_classifier.joblib', 'xgboost_classifier.model']:
        (path / name).write_bytes(b'local model')
    return path


def _fake(monkeypatch: pytest.MonkeyPatch, *, wrong_barcode: bool = False) -> dict:
    calls = {}
    settings = SimpleNamespace(seed=17, cuml=True, accelerator='auto', compute_umap_embedding=True, return_probabilities=False)

    def process(**kwargs):
        calls['process'] = kwargs
        assert kwargs['prediction_mode'] == 'fast'
        assert kwargs['cl_obo_folder'] is False
        # Upstream mutates query/reference and can write auxiliary model files.
        kwargs['query_adata'].obs['mutated'] = True
        kwargs['ref_adata'].obs['mutated'] = True
        (Path(kwargs['save_path_trained_models']) / 'celltypist.pkl').write_bytes(b'changed copy')
        prepared = kwargs['query_adata'][['c1', 'c0']].copy()
        if wrong_barcode:
            prepared.obs_names = ['alien', 'c0']
        prepared.obs['_dataset'] = 'query'
        return SimpleNamespace(adata=prepared)

    def annotate(adata, *, methods):
        assert settings.seed == 0 and settings.cuml is False
        assert settings.accelerator == 'cpu' and settings.compute_umap_embedding is False
        calls['methods'] = methods
        adata.obs['popv_prediction'] = ['B', 'T']
        adata.obs['popv_prediction_score'] = [2, 3]
        adata.obs['popv_celltypist_prediction'] = ['B', 'T']
        adata.obs['popv_svm_prediction'] = ['B', 'T']
        adata.obs['popv_xgboost_prediction'] = ['T', 'T']

    fake = SimpleNamespace(settings=settings, __version__='fake')
    monkeypatch.setitem(sys.modules, 'popv', fake)
    monkeypatch.setitem(sys.modules, 'popv.preprocessing', SimpleNamespace(Process_Query=process))
    monkeypatch.setitem(sys.modules, 'popv.annotation', SimpleNamespace(annotate_data=annotate))
    monkeypatch.setattr(crosschecks.importlib.util, 'find_spec', lambda name: object())
    calls['settings'] = settings
    return calls


def test_popv_preserves_inputs_aligns_predictions_and_retains_filtered_cells(monkeypatch, tmp_path):
    calls = _fake(monkeypatch)
    query = _input()
    reference = _input()
    reference.obs['type'] = ['T', 'B', 'T']
    models = _bundle(tmp_path / 'models')
    cells, status = popv_hypothesis(query, reference, model_dir=models, reference_label='type', cluster_column='cluster', output_dir=tmp_path / 'out')
    assert status['status'] == 'complete'
    assert cells['barcode'].tolist() == ['c0', 'c1', 'c2']
    assert cells['predicted_label'].tolist() == ['T', 'B', 'Unknown']
    assert cells['agreement_count'].tolist() == [3, 2, 0]
    assert cells['abstention_reason'].tolist() == ['', '', 'filtered_by_popv']
    assert cells['popv_celltypist_prediction'].tolist() == ['T', 'B', 'Unknown']
    assert 'mutated' not in query.obs and 'mutated' not in reference.obs
    assert (models / 'celltypist.pkl').read_bytes() == b'local model'
    assert Path(calls['process']['save_path_trained_models']) != models
    assert calls['settings'].seed == 17 and calls['settings'].cuml is True
    assert status['dependent_methods'] == ['celltypist']
    assert status['trained'] is False and status['downloads_requested'] is False
    assert status['filtered_cells'] == 1
    summary = pd.read_csv(status['artifacts']['clusters'], sep='\t')
    assert summary.query("cluster_id == 'a'")['fraction'].tolist() == [0.5, 0.5]
    assert json.loads((tmp_path / 'out/provider_status.json').read_text()) == status


def test_popv_missing_dependency_and_bundle_are_explicit(monkeypatch, tmp_path):
    query = _input()
    reference = _input()
    reference.obs['type'] = ['T', 'B', 'T']
    monkeypatch.setattr(crosschecks.importlib.util, 'find_spec', lambda name: None)
    cells, status = popv_hypothesis(query, reference, model_dir=tmp_path / 'missing', reference_label='type')
    assert cells.empty and status['status'] == 'unavailable'
    assert 'not installed' in status['reason']
    monkeypatch.setattr(crosschecks.importlib.util, 'find_spec', lambda name: object())
    cells, status = popv_hypothesis(query, reference, model_dir=tmp_path / 'missing', reference_label='type')
    assert cells.empty and status['status'] == 'unavailable'
    assert 'bundle misses' in status['reason']


def test_popv_rejects_incompatible_measurements_and_methods(monkeypatch, tmp_path):
    _fake(monkeypatch)
    query, reference = _input(), _input()
    reference.obs['type'] = ['T', 'B', 'T']
    models = _bundle(tmp_path / 'models')
    with pytest.raises(ContractError, match='GPU'):
        popv_hypothesis(query, reference, model_dir=models, reference_label='type', methods=['SCANVI_POPV'])
    query.layers['counts'][0, 0] = 1.5
    with pytest.raises(ContractError, match='integer counts'):
        popv_hypothesis(query, reference, model_dir=models, reference_label='type')
    query = _input()[:, ['G1']].copy()
    _, status = popv_hypothesis(query, reference, model_dir=models, reference_label='type')
    assert status['status'] == 'failed' and 'zero-padding' in status['reason']


def test_popv_never_accepts_unexpected_barcodes(monkeypatch, tmp_path):
    _fake(monkeypatch, wrong_barcode=True)
    reference = _input()
    reference.obs['type'] = ['T', 'B', 'T']
    cells, status = popv_hypothesis(_input(), reference, model_dir=_bundle(tmp_path / 'models'), reference_label='type')
    assert cells.empty and status['status'] == 'failed'
    assert 'unexpected query barcodes' in status['reason']


def test_file_adapter_pins_query_and_reference(monkeypatch, tmp_path):
    _fake(monkeypatch)
    reference = _input()
    reference.obs['type'] = ['T', 'B', 'T']
    _input().write_h5ad(tmp_path / 'query.h5ad')
    reference.write_h5ad(tmp_path / 'ref.h5ad')
    cells, status = predict_popv_file(input_path=tmp_path / 'query.h5ad', reference_path=tmp_path / 'ref.h5ad', model_dir=_bundle(tmp_path / 'models'), reference_label='type', cluster_column='cluster', output_dir=tmp_path / 'out')
    assert len(cells) == 3 and len(status['input_sha256']) == len(status['reference_sha256']) == 64
    assert set(status['artifacts']) == {'cells', 'clusters'}


def test_celltypist_retains_raw_labels_and_realigns_barcodes(monkeypatch, tmp_path):
    model = tmp_path / 'model.pkl'
    model.write_bytes(b'local')
    monkeypatch.setattr(crosschecks.importlib.util, 'find_spec', lambda name: object())
    raw = pd.DataFrame({'predicted_labels': ['B', 'T', 'T'], 'majority_voting': ['B', 'B', 'B']}, index=['c1', 'c2', 'c0'])
    confidence = pd.DataFrame({'conf_score': [.7, .8, .9]}, index=['c1', 'c2', 'c0'])
    result = SimpleNamespace(predicted_labels=raw, to_adata=lambda **kwargs: SimpleNamespace(obs=confidence))
    monkeypatch.setitem(sys.modules, 'celltypist', SimpleNamespace(models=SimpleNamespace(models_path=tmp_path), annotate=lambda *args, **kwargs: result))
    cells, status = celltypist_hypothesis(_input(), model=str(model), cluster_column='cluster')
    assert status['status'] == 'complete'
    assert cells['predicted_label'].tolist() == ['T', 'B', 'T']
    assert cells['majority_voting_label'].tolist() == ['B', 'B', 'B']
    assert cells['uncertainty'].tolist() == pytest.approx([.1, .3, .2])


def test_reconciliation_preserves_marker_calls_and_reports_conflicts():
    markers = pd.DataFrame({'cluster_id': ['a', 'b'], 'label': ['CD8', 'Unknown']})
    membership = pd.DataFrame({'barcode': ['1', '2', '3', '4'], 'cluster_id': ['a', 'a', 'b', 'b']})
    celltypist = membership.assign(predicted_label=['CD8 raw', 'CD8 raw', 'CD8 raw', 'CD8 raw'])
    popv = membership.assign(predicted_label=['B', 'B', 'CD4', 'CD4'], popv_celltypist_prediction='T')
    result = reconcile_predictions(markers, {'popv': popv, 'celltypist': celltypist}, label_mapping={'celltypist': {'CD8 raw': 'CD8'}}, parent_by_label={'CD8': 'T', 'CD4': 'T', 'T': 'Immune', 'B': 'Immune'}, cell_membership=membership)
    assert result['label'].tolist() == markers['label'].tolist()
    assert result['classifier_conflict'].tolist() == [True, False]
    assert result['classifier_conflicting_providers'].tolist() == ['popv', '']
    assert result['classifier_parent_suggestion'].tolist() == ['Immune', 'T']
    assert 'no independence' in result.iloc[0]['classifier_dependency_notes']
    shuffled = reconcile_predictions(markers, {'celltypist': celltypist.iloc[::-1], 'popv': popv.iloc[::-1]}, label_mapping={'celltypist': {'CD8 raw': 'CD8'}}, parent_by_label={'CD8': 'T', 'CD4': 'T', 'T': 'Immune', 'B': 'Immune'}, cell_membership=membership)
    pd.testing.assert_frame_equal(result, shuffled)


def test_reconciliation_rejects_missing_cells_and_keeps_unmapped_labels_unresolved():
    markers = pd.DataFrame({'cluster_id': ['a'], 'label': ['T']})
    cells = pd.DataFrame({'barcode': ['1', '2'], 'cluster_id': ['a', 'a'], 'predicted_label': ['invented', 'invented']})
    result = reconcile_predictions(markers, {'method': cells})
    assert result.iloc[0]['classifier_conflict'] == False  # noqa: E712
    assert json.loads(result.iloc[0]['classifier_summary'])['method']['unresolved'] is True
    with pytest.raises(ContractError, match='barcodes differ'):
        reconcile_predictions(markers, {'method': cells.iloc[:1]}, cell_membership=cells[['barcode', 'cluster_id']])
    with pytest.raises(ContractError, match='cycle'):
        reconcile_predictions(markers, {}, parent_by_label={'T': 'B', 'B': 'T'})
