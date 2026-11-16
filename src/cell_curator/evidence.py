"""Deterministic marker, program, multimodal, and abstention evidence engine."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd
from scipy import optimize, sparse, stats
from scipy.special import expit

from .contracts import (
    ContractError,
    feature_positions,
)
from .models import Capability, MarkerProgram

LOG = logging.getLogger(__name__)

CPU_MARKER_COLUMNS = [
    "cluster_id",
    "method",
    "feature",
    "rank",
    "target_mean",
    "reference_mean",
    "effect",
    "effect_semantics",
    "pvalue",
    "adjusted_pvalue",
    "target_fraction_detected",
    "reference_fraction_detected",
    "method_concordant",
]


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    if values.size == 0:
        return values
    values = np.nan_to_num(values, nan=1.0, posinf=1.0, neginf=1.0)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _matrix_rows(matrix: Any, mask: np.ndarray) -> np.ndarray:
    selected = matrix[mask]
    return selected.toarray() if sparse.issparse(selected) else np.asarray(selected)


def compute_bottom_up_markers(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    *,
    detection: Any | None,
    top_n: int | None = 100,
    marker_config: dict[str, Any] | None = None,
    methods: tuple[str, ...] = ("wilcoxon",),
) -> pd.DataFrame:
    """Compute Wilcoxon, and optionally balanced L2 logistic, one-vs-rest marker evidence.

    Decision row 18: Wilcoxon is the default. The corpus workflow asks for both methods.

    The rank-sum test uses scipy's equivalent Mann-Whitney U test. Logistic
    coefficients come from a multivariate model with an unpenalized intercept;
    they carry no P values. Every observed group gets its own binary fit.
    """

    if matrix.ndim != 2 or matrix.shape[0] != len(groups) or matrix.shape[1] != len(features):
        raise ContractError("matrix, groups, and feature dimensions disagree")
    if (top_n is not None and top_n < 1) or len(features) < 1:
        raise ContractError("marker evidence requires positive top_n and at least one feature")
    if groups.isna().any() or (groups.astype(str).str.strip() == "").any():
        raise ContractError("marker groups must not contain missing or blank labels")
    feature_positions(features)
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all():
        raise ContractError("marker matrix must contain finite values")
    if detection is not None:
        detected = detection.data if sparse.issparse(detection) else np.asarray(detection)
        if detection.shape != matrix.shape or not np.isfinite(detected).all():
            raise ContractError("detection matrix must be finite and match marker dimensions")
        if (detected < 0).any():
            raise ContractError("detection matrix must be nonnegative")
    if "wilcoxon" not in methods or not set(methods) <= {"wilcoxon", "logreg"}:
        raise ContractError("marker methods must include wilcoxon and may add logreg")
    config = marker_config or {}
    wilcoxon = config.get("wilcoxon", {})
    logreg = config.get("logreg", {})
    if logreg.get("class_weight", "balanced") != "balanced":
        raise ContractError("CPU logistic marker evidence requires balanced class weights")
    if logreg.get("penalty", "l2") != "l2":
        raise ContractError("CPU logistic marker evidence requires an L2 penalty")
    regularization = float(logreg.get("C", 1.0))
    tolerance = float(logreg.get("tol", 0.0001))
    max_iter = int(logreg.get("max_iter", 1000))
    retry_max_iter = int(logreg.get("retry_max_iter", 5000))
    if not np.isfinite(regularization) or regularization <= 0:
        raise ContractError("CPU logistic C must be finite and positive")
    if (
        not np.isfinite(tolerance)
        or tolerance <= 0
        or max_iter < 1
        or retry_max_iter < max_iter
    ):
        raise ContractError("CPU logistic tolerance and iteration limits are invalid")
    if wilcoxon.get("correction", "benjamini-hochberg") != "benjamini-hochberg":
        raise ContractError("CPU Wilcoxon requires benjamini-hochberg correction")
    if wilcoxon.get("tie_correct", True) is not True:
        raise ContractError(
            "CPU Mann-Whitney uses tie correction; tie_correct=false is unsupported"
        )
    rows: list[dict[str, Any]] = []
    group_values = groups.astype(str).to_numpy()
    categories = sorted(pd.unique(group_values))
    if len(categories) < 2:
        raise ContractError("CPU marker evidence requires at least two observed groups")
    fit_diagnostics: dict[str, dict[str, Any]] = {}
    for cluster_id in categories:
        target_mask = group_values == cluster_id
        reference_mask = ~target_mask
        if not target_mask.any() or not reference_mask.any():
            continue
        target = _matrix_rows(matrix, target_mask)
        reference = _matrix_rows(matrix, reference_mask)
        target_mean = np.asarray(target.mean(axis=0)).ravel()
        reference_mean = np.asarray(reference.mean(axis=0)).ravel()
        effect = target_mean - reference_mean
        pvalues = []
        for column in range(len(features)):
            try:
                pvalue = stats.mannwhitneyu(
                    target[:, column],
                    reference[:, column],
                    alternative="two-sided",
                    method="asymptotic",
                    use_continuity=bool(wilcoxon.get("continuity_correct", True)),
                ).pvalue
            except ValueError:
                pvalue = 1.0
            pvalues.append(float(pvalue))
        adjusted = benjamini_hochberg(np.asarray(pvalues))
        coefficients = None
        if "logreg" in methods:
            coefficients, diagnostics = _cpu_logistic_markers(
                matrix,
                target_mask,
                C=regularization,
                tolerance=tolerance,
                max_iter=max_iter,
                retry_max_iter=retry_max_iter,
            )
            fit_diagnostics[str(cluster_id)] = diagnostics
        if detection is not None:
            target_detect = (_matrix_rows(detection, target_mask) > 0).mean(axis=0)
            reference_detect = (_matrix_rows(detection, reference_mask) > 0).mean(axis=0)
        else:
            target_detect = np.full(len(features), np.nan)
            reference_detect = np.full(len(features), np.nan)
        for method, scores in (("wilcoxon", effect), ("logreg", coefficients)):
            if scores is None:
                continue
            order = np.lexsort((features.astype(str).to_numpy(), -scores))[:top_n]
            for rank, column in enumerate(order, 1):
                rows.append(
                    {
                        "cluster_id": cluster_id,
                        "method": method,
                        "feature": str(features[column]),
                        "rank": rank,
                        "target_mean": float(target_mean[column]),
                        "reference_mean": float(reference_mean[column]),
                        "effect": float(scores[column]),
                        "effect_semantics": (
                            "mean_difference"
                            if method == "wilcoxon"
                            else "logistic_coefficient"
                        ),
                        "pvalue": float(pvalues[column]) if method == "wilcoxon" else np.nan,
                        "adjusted_pvalue": (
                            float(adjusted[column]) if method == "wilcoxon" else np.nan
                        ),
                        "target_fraction_detected": float(target_detect[column]),
                        "reference_fraction_detected": float(reference_detect[column]),
                    }
                )
    result = pd.DataFrame(rows)
    positive = {
        (str(row.cluster_id), str(row.feature), str(row.method))
        for row in result.itertuples(index=False)
        if float(row.effect) > 0
    }
    result["method_concordant"] = [
        all((str(row.cluster_id), str(row.feature), method) in positive for method in methods)
        for row in result.itertuples(index=False)
    ]
    import scipy

    result.attrs["ranking_group_coverage"] = {method: list(categories) for method in methods}
    result.attrs["runtime"] = {
        "scipy": scipy.__version__,
        "native_api": ["scipy.stats.mannwhitneyu", "scipy.optimize.minimize"],
        "ranking_group_coverage": result.attrs["ranking_group_coverage"],
        **(
            {
                "logreg": {
                    "solver": "L-BFGS-B",
                    "class_weight": "balanced",
                    "penalty": "l2",
                    "C": regularization,
                    "tol": tolerance,
                    "max_iter": max_iter,
                    "retry_max_iter": retry_max_iter,
                    "fit_intercept": True,
                    "objective": "balanced_weighted_sum_log_loss_plus_l2_over_2C",
                    "fits": fit_diagnostics,
                }
            }
            if "logreg" in methods
            else {}
        ),
        "wilcoxon": {
            "tie_correct": True,
            "continuity_correct": bool(wilcoxon.get("continuity_correct", True)),
            "correction": "benjamini-hochberg",
            "alternative": "two-sided",
        },
    }
    return result


def _cpu_logistic_markers(
    matrix: Any,
    target: np.ndarray,
    *,
    C: float,
    tolerance: float,
    max_iter: int,
    retry_max_iter: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Fit the convex balanced logistic objective without penalizing the intercept."""

    labels = target.astype(float)
    sample_weight = np.where(
        target, len(target) / (2 * target.sum()), len(target) / (2 * (~target).sum())
    )

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        coefficients = parameters[:-1]
        linear = np.asarray(matrix @ coefficients).ravel() + parameters[-1]
        residual = sample_weight * (expit(linear) - labels)
        loss = float(np.sum(sample_weight * (np.logaddexp(0, linear) - labels * linear)))
        loss += float(coefficients @ coefficients) / (2 * C)
        gradient = np.r_[
            np.asarray(matrix.T @ residual).ravel() + coefficients / C, residual.sum()
        ]
        return loss, gradient

    fitted = optimize.minimize(
        objective,
        np.zeros(matrix.shape[1] + 1),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": max_iter, "gtol": tolerance, "ftol": np.finfo(float).eps},
    )
    retried = False
    if not fitted.success and retry_max_iter > max_iter:
        retried = True
        fitted = optimize.minimize(
            objective,
            fitted.x,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": retry_max_iter, "gtol": tolerance, "ftol": np.finfo(float).eps},
        )
    if not fitted.success or not np.isfinite(fitted.x).all() or not np.isfinite(fitted.fun):
        raise ContractError(f"CPU logistic regression did not converge: {fitted.message}")
    return np.asarray(fitted.x[:-1]), {
        "converged": True,
        "iterations": int(fitted.nit),
        "retried": retried,
        "intercept": float(fitted.x[-1]),
        "termination": str(fitted.message),
    }


def score_marker_programs(
    matrix: Any,
    groups: pd.Series,
    features: pd.Index,
    programs: list[MarkerProgram],
    *,
    detection: Any | None,
    run_id: str,
    assay: str,
    representation: str,
    capability: Capability,
    positive_detection_threshold: float = 0.1,
    min_negative_coverage: float = 0.35,
) -> pd.DataFrame:
    """Score signed programs without treating unmeasured features as negatives.

    Expected negatives are judged over the ones the panel measures, when at least
    ``min_negative_coverage`` of them are measured; below that the violation is
    unavailable and the gate abstains. Pipelines pass ``min_positive_coverage``, so
    both signs need the same share of their genes measured.

    Genes are standardised across cells, not across cluster means: a z-score over
    the handful of cluster means ignored expression magnitude and moved with the
    cluster composition. Each cell then votes for its best program among the
    siblings sharing a parent and axis, so ``vote_fraction`` is the share of the
    cluster that looks like the program, which is the cluster's purity for that
    label. The expected-negative test is run on the cells that voted for the
    program: a neighbouring lineage's cells in an impure cluster no longer veto
    the call, while a negative marker expressed by the supporting cells still does.
    """

    group_values = groups.astype(str).to_numpy()
    cluster_ids = sorted(set(group_values))
    cluster_index = {cluster: index for index, cluster in enumerate(cluster_ids)}
    cell_cluster = np.asarray([cluster_index[value] for value in group_values])
    n_clusters = len(cluster_ids)
    cluster_sizes = np.bincount(cell_cluster, minlength=n_clusters)
    positions = feature_positions(features)
    compatible_spaces = {capability.feature_space}
    if capability.assay_kind.value == "gene_activity":
        compatible_spaces.add("gene")
    scored: list[
        tuple[MarkerProgram, list[int], list[int], tuple[str, ...], tuple[str, ...]]
    ] = []
    for program in programs:
        if program.feature_space not in compatible_spaces:
            continue
        positive = [positions[item] for item in program.positive if item in positions]
        negative = [positions[item] for item in program.negative if item in positions]
        if not positive:
            continue
        scored.append(
            (
                program,
                positive,
                negative,
                tuple(item for item in program.positive if item not in positions),
                tuple(item for item in program.negative if item not in positions),
            )
        )
    columns_used = sorted({index for _, pos, neg, _, _ in scored for index in (*pos, *neg)})
    local = {feature: offset for offset, feature in enumerate(columns_used)}
    n_cells = len(group_values)
    width = len(columns_used)
    cell_sum = np.zeros(width)
    cell_square = np.zeros(width)
    cluster_sum = np.zeros((n_clusters, width))
    cluster_detected = np.zeros((n_clusters, width))
    chunk = 20_000

    def block(source: Any, start: int, stop: int) -> np.ndarray:
        rows = source[start:stop]
        rows = rows[:, columns_used] if width else rows[:, :0]
        return rows.toarray() if sparse.issparse(rows) else np.asarray(rows, dtype=float)

    for start in range(0, n_cells, chunk):
        stop = min(start + chunk, n_cells)
        values = block(matrix, start, stop).astype(float)
        cell_sum += values.sum(axis=0)
        cell_square += (values**2).sum(axis=0)
        np.add.at(cluster_sum, cell_cluster[start:stop], values)
        if detection is not None:
            np.add.at(
                cluster_detected, cell_cluster[start:stop], block(detection, start, stop) > 0
            )
    mean = cell_sum / max(n_cells, 1)
    spread = np.sqrt(np.maximum(cell_square / max(n_cells, 1) - mean**2, 0.0))
    spread = np.where(spread <= 1e-12, 1.0, spread)
    sizes = np.maximum(cluster_sizes, 1)[:, None]
    standardized = (cluster_sum / sizes - mean) / spread
    fractions = cluster_detected / sizes if detection is not None else None

    # Per-cell votes among siblings: programs compete only with those that share
    # their parent and axis, never with their own parent or with cell states.
    competitions: dict[tuple[str, str], list[int]] = {}
    for index, (program, *_rest) in enumerate(scored):
        competitions.setdefault((program.parent or "", program.axis), []).append(index)
    votes = np.zeros((n_clusters, len(scored)))
    voter_detected = [np.zeros((n_clusters, len(neg))) for _, _, neg, _, _ in scored]
    for start in range(0, n_cells, chunk):
        stop = min(start + chunk, n_cells)
        z = (block(matrix, start, stop) - mean) / spread
        detected = block(detection, start, stop) > 0 if detection is not None else None
        clusters = cell_cluster[start:stop]
        for members in competitions.values():
            cell_scores = np.column_stack(
                [
                    z[:, [local[i] for i in scored[m][1]]].mean(axis=1)
                    - (
                        np.maximum(z[:, [local[i] for i in scored[m][2]]], 0.0).mean(axis=1)
                        if scored[m][2]
                        else 0.0
                    )
                    for m in members
                ]
            )
            winner = np.asarray(members)[np.argmax(cell_scores, axis=1)]
            # A cell no program fits better than the dataset average does not vote,
            # but stays in the cluster size, so a low-signal cluster loses support
            # instead of handing its cells to whichever sibling sorts first.
            winner[cell_scores.max(axis=1) <= 0] = -1
            for m in members:
                voted = winner == m
                np.add.at(votes[:, m], clusters[voted], 1)
                if detected is not None and scored[m][2]:
                    np.add.at(
                        voter_detected[m],
                        clusters[voted],
                        detected[np.ix_(voted, [local[i] for i in scored[m][2]])],
                    )

    rows: list[dict[str, Any]] = []
    for m, (program, positive, negative, missing_positive, missing_negative) in enumerate(
        scored
    ):
        program_coverage = len(positive) / len(program.positive)
        negative_coverage = len(negative) / len(program.negative) if program.negative else 1.0
        pos_local = [local[i] for i in positive]
        neg_local = [local[i] for i in negative]
        for row_index, cluster_id in enumerate(cluster_ids):
            positive_signal = float(np.mean(standardized[row_index, pos_local]))
            negative_signal = (
                float(np.mean(np.maximum(standardized[row_index, neg_local], 0.0)))
                if negative
                else 0.0
            )
            voters = votes[row_index, m]
            if fractions is not None:
                positive_coverage = float(
                    np.mean(fractions[row_index, pos_local] >= positive_detection_threshold)
                )
                if negative_coverage < min_negative_coverage:
                    negative_violation = float("nan")
                elif negative:
                    observed = (
                        voter_detected[m][row_index] / voters
                        if voters > 0
                        else fractions[row_index, neg_local]
                    )
                    negative_violation = float(
                        np.mean(observed >= positive_detection_threshold)
                    )
                else:
                    negative_violation = 0.0
            else:
                positive_coverage = float("nan")
                negative_violation = float("nan")
            score = (
                positive_signal
                - negative_signal
                - (negative_violation if np.isfinite(negative_violation) else 0.0)
            )
            rows.append(
                {
                    "schema_version": 3,
                    "run_id": run_id,
                    "cluster_id": cluster_id,
                    "label": program.label,
                    "axis": program.axis,
                    "parent_label": program.parent or "",
                    "cl_id": program.cl_id or "",
                    "assay": assay,
                    "representation": representation,
                    "capability_level": capability.level.value,
                    "feature_space": capability.feature_space,
                    "score": score,
                    "positive_signal": positive_signal,
                    "negative_signal": negative_signal,
                    "positive_coverage": positive_coverage,
                    "negative_violation": negative_violation,
                    "program_coverage": program_coverage,
                    "measured_positive": len(positive),
                    "total_positive": len(program.positive),
                    "measured_negative": len(negative),
                    "total_negative": len(program.negative),
                    "missing_positive": ";".join(missing_positive),
                    "missing_negative": ";".join(missing_negative),
                    "detection_valid": detection is not None,
                    "source": program.source,
                    "vote_fraction": float(voters / max(cluster_sizes[row_index], 1)),
                    "n_voters": int(voters),
                }
            )
    columns = [
        "schema_version",
        "run_id",
        "cluster_id",
        "label",
        "axis",
        "parent_label",
        "cl_id",
        "assay",
        "representation",
        "capability_level",
        "feature_space",
        "score",
        "positive_signal",
        "negative_signal",
        "positive_coverage",
        "negative_violation",
        "program_coverage",
        "measured_positive",
        "total_positive",
        "measured_negative",
        "total_negative",
        "missing_positive",
        "missing_negative",
        "detection_valid",
        "source",
        "vote_fraction",
        "n_voters",
    ]
    return pd.DataFrame(rows, columns=columns)
