#!/usr/bin/env python3
"""Analyze structural spectral follow-ups with matched controls and diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.multitest import multipletests


SEED = 20260922
TARGET_LAMBDA2 = 2 + 2 * np.sqrt(2)
COLORS = {5: "#0072B2", 6: "#E69F00", 7: "#009E73", 8: "#CC79A7"}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=True)


def bootstrap_mean_ci(values: np.ndarray, rng: np.random.Generator, reps: int = 5000) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    estimates = np.empty(reps)
    for start in range(0, reps, 250):
        count = min(250, reps - start)
        estimates[start:start + count] = rng.choice(values, (count, len(values)), replace=True).mean(axis=1)
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def prepare_cases(census: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    cases = read_tsv(census / "SPECTRUM.tsv")
    cases = cases[cases["kind"] == "exact_component"].copy()
    cases["group_degree"] = cases["group"].astype(str) + cases["degree"].astype(str)
    cases["relaxation_gap"] = cases["normalized_gap"]
    cases["absolute_gap"] = 1 - cases["abs_nontrivial_radius"] / 6
    cases["log_vertices"] = np.log(cases["quotient_states"].clip(lower=1))
    cases["loop_rate"] = cases["loop_moves"] / (6 * cases["quotient_states"])
    cases["parallel_rate"] = cases["parallel_excess"] / (6 * cases["quotient_states"])
    # S8 comes from a convenience sample of completed searches, unlike the
    # stored exact component catalogs for the other strata.
    cases["sampling_frame"] = np.where((cases["group"] == "S") & (cases["degree"] == 8), "sampled_S8", "stored_exact_catalog")

    diagnostic_path = census / "GRAPH_DIAGNOSTICS.tsv"
    if diagnostic_path.exists():
        diagnostics = read_tsv(diagnostic_path)
        actual = diagnostics[diagnostics["source"] == "case"].copy()
        keep = ["case_id", "generator_1_cycles", "generator_1_fixed", "generator_2_fixed", "generator_3_fixed",
                "generator_1_cycle_count", "generator_2_cycle_count", "generator_3_cycle_count",
                "braid_12_defects", "braid_23_defects", "commute_13_defects",
                "exact_bipartite", "maxcut_conflict_fraction", "wl_hash"]
        cases = cases.merge(actual[keep], on="case_id", how="left")
        cases["generator_cycle_density"] = cases["generator_1_cycle_count"] / cases["quotient_states"]
    else:
        diagnostics = pd.DataFrame()

    baselines = read_tsv(census / "BASELINES.tsv")
    baselines["relaxation_gap"] = baselines["normalized_gap"]
    baselines["absolute_gap"] = 1 - baselines["abs_nontrivial_radius"] / 6
    means = baselines.groupby(["case_id", "baseline"])[["relaxation_gap", "absolute_gap"]].mean().reset_index()
    for outcome in ("relaxation_gap", "absolute_gap"):
        wide = means.pivot(index="case_id", columns="baseline", values=outcome)
        for kind in wide.columns:
            cases = cases.merge(wide[[kind]].rename(columns={kind: f"baseline_{outcome}_{kind}"}), left_on="case_id", right_index=True, how="left")
            cases[f"excess_{outcome}_{kind}"] = cases[outcome] - cases[f"baseline_{outcome}_{kind}"]
    return cases, diagnostics


def design_matrix(frame: pd.DataFrame, include_genus: bool) -> pd.DataFrame:
    categorical = pd.get_dummies(frame["group_degree"], prefix="stratum", dtype=float, drop_first=True)
    numeric_names = ["log_vertices", "loop_rate", "parallel_rate"]
    # Braid generators are conjugate in B4, so all three fixed-point and
    # cycle-count columns coincide on genuine Hurwitz actions. Moreover,
    # fixed points / n equal the total loop rate exactly. Including those
    # redundant columns created the singular regression in the 2026-09-22 run.
    if "generator_cycle_density" in frame and frame["generator_cycle_density"].notna().all():
        numeric_names.append("generator_cycle_density")
    numeric = frame[numeric_names].astype(float).copy()
    scale = numeric.std(ddof=0).replace(0, 1)
    numeric = (numeric - numeric.mean()) / scale
    if include_genus:
        numeric["genus"] = frame["genus"].astype(float) - frame["genus"].mean()
    matrix = sm.add_constant(pd.concat([categorical.reset_index(drop=True), numeric.reset_index(drop=True)], axis=1), has_constant="add")
    values = matrix.to_numpy(float)
    rank = np.linalg.matrix_rank(values)
    condition = np.linalg.cond(values)
    if rank != values.shape[1] or not np.isfinite(condition) or condition > 1e9:
        raise ValueError(f"nonidentifiable design: rank {rank}/{values.shape[1]}, condition {condition:.3g}, columns={list(matrix)}")
    return matrix


def genus_models(cases: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cases = cases[cases["sampling_frame"] == "stored_exact_catalog"].copy()
    outcomes = ["relaxation_gap", "absolute_gap"]
    outcomes += sorted(name for name in cases if name.startswith("excess_relaxation_gap_"))
    outcomes += sorted(name for name in cases if name.startswith("excess_absolute_gap_"))
    rows, diagnostics = [], []
    for outcome in outcomes:
        frame = cases.dropna(subset=[outcome]).reset_index(drop=True)
        base_x = design_matrix(frame, False)
        full_x = design_matrix(frame, True)
        y = frame[outcome].astype(float)
        cluster = frame["group_degree"].astype(str) + "|" + frame["passport"].astype(str)
        base = sm.OLS(y, base_x).fit()
        model = sm.OLS(y, full_x).fit(cov_type="cluster", cov_kwds={"groups": cluster})
        lo, hi = model.conf_int().loc["genus"]
        rows.append({
            "outcome": outcome, "n": len(frame), "genus_coefficient": model.params["genus"],
            "passport_cluster_se": model.bse["genus"], "ci_low": lo, "ci_high": hi,
            "p_value": model.pvalues["genus"], "r_squared": model.rsquared,
            "adjusted_r_squared": model.rsquared_adj, "delta_r_squared": model.rsquared - base.rsquared,
            "passport_clusters": cluster.nunique(),
        })
        jb = stats.jarque_bera(model.resid)
        bp = het_breuschpagan(model.resid, np.asarray(full_x, dtype=float))
        diagnostics.append({
            "outcome": outcome, "residual_skew": stats.skew(model.resid),
            "jarque_bera": jb.statistic, "jarque_bera_p": jb.pvalue,
            "breusch_pagan_lm": bp[0], "breusch_pagan_p": bp[1],
            "condition_number": np.linalg.cond(np.asarray(full_x, dtype=float)),
            "design_rank": np.linalg.matrix_rank(np.asarray(full_x, dtype=float)),
            "design_columns": full_x.shape[1],
            "inference": "OLS coefficient with group-degree/passport cluster-robust standard error; exploratory observational association",
        })
    result = pd.DataFrame(rows)
    if len(result):
        result["holm_p"] = multipletests(result["p_value"], method="holm")[1]
    return result, pd.DataFrame(diagnostics)


def passport_cross_validation(cases: pd.DataFrame) -> pd.DataFrame:
    cases = cases[cases["sampling_frame"] == "stored_exact_catalog"].copy()
    rows = []
    for outcome in ("relaxation_gap", "absolute_gap"):
        frame = cases.dropna(subset=[outcome]).reset_index(drop=True)
        matrices = {False: np.asarray(design_matrix(frame, False), dtype=float), True: np.asarray(design_matrix(frame, True), dtype=float)}
        y = frame[outcome].to_numpy(float)
        passports = frame["passport"].astype(str).to_numpy()
        unique = np.unique(passports)
        predictions = {key: np.full(len(frame), np.nan) for key in matrices}
        for passport in unique:
            test = passports == passport
            train = ~test
            if any(train.sum() <= x.shape[1] or np.linalg.matrix_rank(x[train]) < x.shape[1] for x in matrices.values()):
                continue
            for key, x in matrices.items():
                beta = np.linalg.lstsq(x[train], y[train], rcond=None)[0]
                predictions[key][test] = x[test] @ beta
        valid = np.isfinite(predictions[False]) & np.isfinite(predictions[True])
        if not valid.any():
            raise ValueError(f"no identifiable held-out passports for {outcome}")
        for include_genus in matrices:
            prediction = predictions[include_genus]
            errors = prediction[valid] - y[valid]
            rows.append({
                "outcome": outcome, "model": "structural_plus_genus" if include_genus else "structural",
                "n": valid.sum(), "passports": len(unique), "rmse": np.sqrt(np.mean(errors ** 2)),
                "mae": np.mean(np.abs(errors)), "r_squared_prediction": 1 - np.sum(errors ** 2) / np.sum((y[valid] - y[valid].mean()) ** 2),
            })
    result = pd.DataFrame(rows)
    for outcome in result["outcome"].unique():
        block = result[result["outcome"] == outcome].set_index("model")
        delta = block.loc["structural_plus_genus", "rmse"] - block.loc["structural", "rmse"]
        result.loc[result["outcome"] == outcome, "genus_delta_rmse"] = delta
    return result


def matched_hurwitz_slopes(cases: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cases = cases[cases["sampling_frame"] == "stored_exact_catalog"].copy()
    match_features = ["log_vertices", "loop_rate", "parallel_rate"]
    match_features += [name for name in ("generator_cycle_density",) if name in cases and cases[name].notna().all()]
    exact_profile = "generator_1_cycles" in cases and cases["generator_1_cycles"].notna().all()
    records = []
    for _, block in cases.groupby("group_degree"):
        block = block.reset_index(drop=True)
        values = block[match_features].to_numpy(float)
        values = (values - values.mean(axis=0)) / np.where(values.std(axis=0) == 0, 1, values.std(axis=0))
        for i in range(len(block)):
            candidate = block["genus"].to_numpy() != block.loc[i, "genus"]
            tier = "group_degree"
            if exact_profile:
                exact = candidate & (block["generator_1_cycles"].to_numpy() == block.loc[i, "generator_1_cycles"])
                if exact.any():
                    candidate = exact
                    tier = "generator_cycle_type"
            indices = np.flatnonzero(candidate)
            if not len(indices):
                continue
            distance = ((values[indices] - values[i]) ** 2).sum(axis=1)
            j = int(indices[np.argmin(distance)])
            first, second = (i, j) if block.loc[i, "genus"] < block.loc[j, "genus"] else (j, i)
            records.append({
                "case_low": block.loc[first, "case_id"], "case_high": block.loc[second, "case_id"],
                "group_degree": block.loc[first, "group_degree"], "genus_low": block.loc[first, "genus"],
                "genus_high": block.loc[second, "genus"], "match_distance": float(np.sqrt(distance.min())), "match_tier": tier,
                "relaxation_slope": (block.loc[second, "relaxation_gap"] - block.loc[first, "relaxation_gap"]) / (block.loc[second, "genus"] - block.loc[first, "genus"]),
                "absolute_slope": (block.loc[second, "absolute_gap"] - block.loc[first, "absolute_gap"]) / (block.loc[second, "genus"] - block.loc[first, "genus"]),
            })
    # A case can have many nearest-neighbor suitors. Keep only disjoint pairs;
    # their summaries remain descriptive because passports repeat graph types.
    used = set()
    disjoint = []
    for row in sorted(records, key=lambda x: (x["match_tier"] != "generator_cycle_type", x["match_distance"])):
        if row["case_low"] not in used and row["case_high"] not in used:
            used.update((row["case_low"], row["case_high"]))
            disjoint.append(row)
    pair_frame = pd.DataFrame(disjoint)
    if pair_frame.empty:
        return pair_frame, pd.DataFrame()
    summaries = []
    for outcome in ("relaxation_slope", "absolute_slope"):
        values = pair_frame[outcome].to_numpy(float)
        summaries.append({"outcome": outcome, "pairs": len(values), "mean_slope_per_genus": values.mean(), "median_slope_per_genus": np.median(values), "exact_profile_fraction": np.mean(pair_frame["match_tier"] == "generator_cycle_type"), "inference": "descriptive only; disjoint cases but recurring passports and graph types"})
    return pair_frame, pd.DataFrame(summaries)


def relation_summary(diagnostics: pd.DataFrame) -> pd.DataFrame:
    if diagnostics.empty:
        return pd.DataFrame()
    rows = []
    for source, block in diagnostics.groupby("source"):
        any_defect = (block[["braid_12_defects", "braid_23_defects", "commute_13_defects"]].sum(axis=1) > 0)
        rows.append({
            "source": source, "graphs": len(block), "fraction_with_relation_defect": any_defect.mean(),
            "mean_braid_12_defect_fraction": (block["braid_12_defects"] / block["vertices"]).mean(),
            "mean_braid_23_defect_fraction": (block["braid_23_defects"] / block["vertices"]).mean(),
            "mean_commute_13_defect_fraction": (block["commute_13_defects"] / block["vertices"]).mean(),
        })
    return pd.DataFrame(rows)


def bipartite_bounds(cases: pd.DataFrame, deep: Path | None = None) -> pd.DataFrame:
    """Bound minimum frustrated directed-edge density without calling it a gap."""
    columns = ["case_id", "group_degree", "genus", "sampling_frame", "quotient_states",
               "loop_rate", "absolute_gap", "method", "lambda_min"]
    if "exact_bipartite" not in cases:
        return pd.DataFrame()
    columns += ["exact_bipartite"]
    frame = cases[columns].copy()
    frame["frustration_upper_bound_heuristic"] = np.nan
    if deep is not None and (deep / "GRAPH_DIAGNOSTICS.tsv").exists():
        diagnostics = read_tsv(deep / "GRAPH_DIAGNOSTICS.tsv")
        actual = diagnostics[diagnostics["source"] == "case"][["case_id", "maxcut_conflict_fraction"]]
        frame = frame.merge(actual.rename(columns={"maxcut_conflict_fraction": "frustration_upper_bound_heuristic"}), on="case_id", how="left", validate="one_to_one", suffixes=("", "_deep"))
        if "frustration_upper_bound_heuristic_deep" in frame:
            frame["frustration_upper_bound_heuristic"] = frame["frustration_upper_bound_heuristic_deep"]
            frame = frame.drop(columns="frustration_upper_bound_heuristic_deep")
    # Jacobi eigenvalues are floating-point, so only the loop term is a
    # rigorous bound here; the spectral term is a useful numerical estimate.
    frame["spectral_lower_bound_estimate"] = np.where(frame["method"] == "jacobi-exact", (1 + frame["lambda_min"] / 6) / 2, np.nan)
    frame["frustration_lower_bound"] = frame["loop_rate"]
    frame["maxcut_evaluated"] = frame["frustration_upper_bound_heuristic"].notna()
    frame["bound_order_valid"] = pd.Series(pd.NA, index=frame.index, dtype="boolean")
    evaluated = frame["maxcut_evaluated"]
    frame.loc[evaluated, "bound_order_valid"] = frame.loc[evaluated, "frustration_lower_bound"] <= frame.loc[evaluated, "frustration_upper_bound_heuristic"] + 1e-8
    return frame


def trace_analysis(deep: Path | None, rng: np.random.Generator) -> pd.DataFrame:
    if deep is None or not (deep / "TRACE_MOMENTS.tsv").exists():
        return pd.DataFrame()
    moments = read_tsv(deep / "TRACE_MOMENTS.tsv")
    actual = moments[moments["source"] == "case"][["case_id", "order", "normalized_trace"]].rename(columns={"normalized_trace": "actual"})
    baseline = moments[moments["source"] != "case"].groupby(["case_id", "source", "order"])["normalized_trace"].mean().reset_index(name="baseline")
    merged = baseline.merge(actual, on=["case_id", "order"])
    metadata = read_tsv(deep / "SPECTRUM.tsv")[["case_id", "group", "degree", "passport"]]
    merged = merged.merge(metadata, on="case_id", validate="many_to_one")
    merged = merged[(merged["group"] != "S") | (merged["degree"] != 8)].copy()
    merged["passport_cluster"] = merged["group"].astype(str) + merged["degree"].astype(str) + "|" + merged["passport"].astype(str)
    merged["difference"] = merged["actual"] - merged["baseline"]
    rows = []
    for (source, order), block in merged.groupby(["source", "order"]):
        # A passport may contribute many related graphs, so use one paired
        # difference per passport cluster for uncertainty and the sign test.
        cluster_means = block.groupby("passport_cluster")[["actual", "baseline", "difference"]].mean()
        difference = cluster_means["difference"].to_numpy(float)
        lo, hi = bootstrap_mean_ci(difference, rng)
        # Zeros carry no direction. The former >=0 / <=0 count counted each
        # zero twice and produced incorrect sign-test p-values.
        nonzero = difference[np.abs(difference) > 1e-10]
        positive = np.count_nonzero(nonzero > 0)
        negative = np.count_nonzero(nonzero < 0)
        p = stats.binomtest(positive, positive + negative, 0.5).pvalue if len(nonzero) else 1.0
        rows.append({
            "baseline": source, "order": order, "n_cases": len(block), "n_passport_clusters": len(difference), "actual_mean": cluster_means["actual"].mean(),
            "baseline_mean": cluster_means["baseline"].mean(), "paired_difference": difference.mean(),
            "ci_low": lo, "ci_high": hi, "standardized_paired_effect": difference.mean() / difference.std(ddof=1) if len(difference) > 1 and difference.std(ddof=1) else np.nan,
            "positive_pairs": positive, "negative_pairs": negative, "ties": len(difference) - len(nonzero),
            "sign_p": p,
        })
    result = pd.DataFrame(rows)
    if len(result):
        result["holm_p"] = multipletests(result["sign_p"], method="holm")[1]
    return result


def convergence_analysis(deep: Path | None) -> pd.DataFrame:
    if deep is None or not (deep / "LANCZOS_CONVERGENCE.tsv").exists():
        return pd.DataFrame()
    frame = read_tsv(deep / "LANCZOS_CONVERGENCE.tsv")
    if frame.empty:
        return frame
    reference = frame.sort_values("steps").groupby("case_id").tail(1)[["case_id", "lambda2", "lambda_min", "steps"]].rename(columns={"lambda2": "reference_lambda2", "lambda_min": "reference_lambda_min", "steps": "reference_steps"})
    frame = frame.merge(reference, on="case_id")
    frame["lambda2_error"] = (frame["lambda2"] - frame["reference_lambda2"]).abs()
    frame["lambda_min_error"] = (frame["lambda_min"] - frame["reference_lambda_min"]).abs()
    return frame.groupby(["steps", "reference_steps"]).agg(cases=("case_id", "size"), median_lambda2_error=("lambda2_error", "median"), max_lambda2_error=("lambda2_error", "max"), median_lambda_min_error=("lambda_min_error", "median"), max_lambda_min_error=("lambda_min_error", "max")).reset_index()


def high_replicate_analysis(path: Path | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    if path is None or not (path / "BASELINES.tsv").exists():
        return pd.DataFrame(), pd.DataFrame()
    spectrum = read_tsv(path / "SPECTRUM.tsv").set_index("case_id")
    baseline = read_tsv(path / "BASELINES.tsv")
    baseline["relaxation_gap"] = baseline["normalized_gap"]
    baseline["absolute_gap"] = 1 - baseline["abs_nontrivial_radius"] / 6
    rows = []
    for (case_id, kind), block in baseline.groupby(["case_id", "baseline"]):
        if case_id not in spectrum.index:
            continue
        actual_relaxation = float(spectrum.loc[case_id, "normalized_gap"])
        actual_absolute = 1 - float(spectrum.loc[case_id, "abs_nontrivial_radius"]) / 6
        for outcome, actual in (("relaxation_gap", actual_relaxation), ("absolute_gap", actual_absolute)):
            values = block[outcome].dropna().to_numpy(float)
            sd = values.std(ddof=1) if len(values) > 1 else np.nan
            rows.append({
                "case_id": case_id, "baseline": kind, "outcome": outcome, "replicates": len(values),
                "actual": actual, "baseline_mean": values.mean(), "actual_minus_baseline": actual - values.mean(),
                "standardized_difference": (actual - values.mean()) / sd if sd > 0 else np.nan,
                "one_sided_empirical_p_actual_smaller": (np.count_nonzero(values <= actual) + 1) / (len(values) + 1),
            })
    detail = pd.DataFrame(rows)
    if detail.empty:
        return detail, pd.DataFrame()
    summary = detail.groupby(["baseline", "outcome"]).agg(cases=("case_id", "size"), mean_difference=("actual_minus_baseline", "mean"), median_standardized_difference=("standardized_difference", "median"), fraction_empirical_p_le_005=("one_sided_empirical_p_actual_smaller", lambda x: np.mean(x <= .05))).reset_index()
    return detail, summary


def switch_sensitivity(specifications: list[str]) -> pd.DataFrame:
    frames = []
    for specification in specifications:
        label, raw_path = specification.split("=", 1)
        path = Path(raw_path).resolve()
        baseline = read_tsv(path / "BASELINES.tsv")
        baseline = baseline[baseline["baseline"] == "random_6_regular"].copy()
        baseline["absolute_gap"] = 1 - baseline["abs_nontrivial_radius"] / 6
        means = baseline.groupby("case_id")[["normalized_gap", "absolute_gap"]].mean().reset_index()
        means["setting"] = label
        frames.append(means)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    reference_label = frames[-1]["setting"].iloc[0]
    reference = combined[combined["setting"] == reference_label].set_index("case_id")
    rows = []
    for setting, block in combined.groupby("setting"):
        joined = block.set_index("case_id").join(reference[["normalized_gap", "absolute_gap"]], rsuffix="_reference", how="inner")
        rows.append({
            "setting": setting, "reference": reference_label, "cases": len(joined),
            "mean_relaxation_gap": joined["normalized_gap"].mean(),
            "mean_difference_from_reference": (joined["normalized_gap"] - joined["normalized_gap_reference"]).mean(),
            "max_absolute_difference_from_reference": (joined["normalized_gap"] - joined["normalized_gap_reference"]).abs().max(),
            "mean_absolute_gap": joined["absolute_gap"].mean(),
            "mean_absolute_gap_difference_from_reference": (joined["absolute_gap"] - joined["absolute_gap_reference"]).mean(),
        })
    return pd.DataFrame(rows)


def recurring_eigenvalue(cases: pd.DataFrame, census: Path) -> pd.DataFrame:
    recurring = cases[(cases["method"] == "jacobi-exact") & ((cases["lambda2"] - TARGET_LAMBDA2).abs() < 1e-8)].copy()
    if recurring.empty:
        return recurring
    if "wl_hash" not in recurring:
        recurring["wl_hash"] = "unavailable"
    rows = []
    values_path = census / "EIGENVALUES.tsv"
    spectral_hashes: dict[str, str] = {}
    if values_path.exists():
        values = read_tsv(values_path)
        values = values[(values["source"] == "case") & values["case_id"].isin(recurring["case_id"])]
        for case_id, block in values.groupby("case_id"):
            packed = ",".join(f"{value:.8f}" for value in block.sort_values("eigen_index")["eigenvalue"])
            spectral_hashes[str(case_id)] = hashlib.sha256(packed.encode()).hexdigest()[:16]
    recurring["spectral_hash"] = recurring["case_id"].map(spectral_hashes).fillna("unavailable")
    for keys, block in recurring.groupby(["quotient_states", "wl_hash", "spectral_hash"], dropna=False):
        rows.append({"vertices": keys[0], "wl_hash": keys[1], "spectral_hash": keys[2], "cases": len(block), "groups": ",".join(sorted(block["group_degree"].unique())), "genera": ",".join(map(str, sorted(block["genus"].unique()))), "passports": block["passport"].nunique()})
    return pd.DataFrame(rows).sort_values("cases", ascending=False)


def make_figure(output: Path, cases: pd.DataFrame, models: pd.DataFrame, traces: pd.DataFrame, convergence: pd.DataFrame, recurring: pd.DataFrame) -> None:
    mpl.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 6.2), layout="constrained")
    rng = np.random.default_rng(SEED)
    for degree, block in cases.groupby("degree"):
        sample = block if len(block) <= 1200 else block.iloc[rng.choice(len(block), 1200, replace=False)]
        axes[0, 0].scatter(sample["relaxation_gap"], sample["absolute_gap"], s=6, alpha=.25, color=COLORS.get(int(degree), "#777777"), label=f"d={degree}")
    axes[0, 0].plot([0, .4], [0, .4], ls="--", lw=.8, color="#333333")
    axes[0, 0].set(xlabel="Relaxation gap", ylabel="Absolute gap", title="A  Negative spectrum controls mixing")
    axes[0, 0].legend(frameon=False, ncol=2)

    rate = cases.assign(near=cases["absolute_gap"] < .01).groupby("group_degree")["near"].mean()
    axes[0, 1].bar(rate.index, rate.values, color=[COLORS.get(int(label[1:]), "#777777") for label in rate.index])
    axes[0, 1].tick_params(axis="x", rotation=45)
    axes[0, 1].set(ylabel="Fraction with absolute gap < 0.01", title="B  Near-bipartite components")

    show = models[models["outcome"].isin(["relaxation_gap", "absolute_gap"])].copy()
    if len(show):
        y = np.arange(len(show))
        axes[0, 2].errorbar(show["genus_coefficient"], y, xerr=[show["genus_coefficient"] - show["ci_low"], show["ci_high"] - show["genus_coefficient"]], fmt="o", color="#0072B2", capsize=3)
        axes[0, 2].axvline(0, color="#333333", lw=.8, ls="--")
        axes[0, 2].set(yticks=y, yticklabels=show["outcome"].str.replace("_", " "), xlabel="Adjusted change per genus", title="C  Genus after structural adjustment")

    if len(traces):
        for source, block in traces.groupby("baseline"):
            axes[1, 0].plot(block["order"], block["paired_difference"], marker="o", ms=3, label=source.replace("_", " "))
        axes[1, 0].axhline(0, color="#333333", lw=.8)
        axes[1, 0].set(xlabel="Moment order k", ylabel="Hurwitz − null normalized trace", title="D  Closed-walk moment excess")
        axes[1, 0].legend(frameon=False)
    else:
        axes[1, 0].text(.5, .5, "Deep run pending", ha="center", va="center", transform=axes[1, 0].transAxes)
        axes[1, 0].set(title="D  Closed-walk moment excess")

    if len(convergence):
        axes[1, 1].semilogy(convergence["steps"], convergence["max_lambda2_error"].clip(lower=1e-12), marker="o", label=r"$\lambda_2$")
        axes[1, 1].semilogy(convergence["steps"], convergence["max_lambda_min_error"].clip(lower=1e-12), marker="s", label=r"$\lambda_{min}$")
        axes[1, 1].set(xlabel="Lanczos steps", ylabel="Maximum error vs longest run", title="E  Numerical convergence")
        axes[1, 1].legend(frameon=False)
    else:
        axes[1, 1].text(.5, .5, "Convergence run pending", ha="center", va="center", transform=axes[1, 1].transAxes)
        axes[1, 1].set(title="E  Numerical convergence")

    if len(recurring):
        by_size = recurring.groupby("vertices")["cases"].sum().sort_values(ascending=False).head(10)
        axes[1, 2].bar([str(x) for x in by_size.index], by_size.values, color="#CC79A7")
        axes[1, 2].tick_params(axis="x", rotation=45)
        axes[1, 2].set(xlabel="Vertices", ylabel="Cases", title=r"F  Near $\lambda_2=2+2\sqrt{2}$ families")
    fig.suptitle("Hurwitz four-tuple spectral follow-up tests", fontsize=12)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output / f"spectral-followups.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--deep", type=Path)
    parser.add_argument("--convergence", type=Path)
    parser.add_argument("--high-reps", type=Path)
    parser.add_argument("--switch-run", action="append", default=[], help="LABEL=PATH; repeat for switch-chain settings")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    census = args.census.resolve(); deep = args.deep.resolve() if args.deep else None
    convergence_root = args.convergence.resolve() if args.convergence else deep
    high_reps_root = args.high_reps.resolve() if args.high_reps else deep
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    cases, diagnostics = prepare_cases(census)
    models, model_diagnostics = genus_models(cases)
    cross_validation = passport_cross_validation(cases)
    pairs, pair_summary = matched_hurwitz_slopes(cases)
    relations = relation_summary(diagnostics)
    bounds = bipartite_bounds(cases, deep)
    traces = trace_analysis(deep, rng)
    convergence = convergence_analysis(convergence_root)
    high_reps, high_reps_summary = high_replicate_analysis(high_reps_root)
    switches = switch_sensitivity(args.switch_run)
    recurring = recurring_eigenvalue(cases, census)

    cases.groupby(["group", "degree", "genus"]).agg(cases=("case_id", "size"), median_vertices=("quotient_states", "median"), mean_relaxation_gap=("relaxation_gap", "mean"), median_relaxation_gap=("relaxation_gap", "median"), mean_absolute_gap=("absolute_gap", "mean"), near_bipartite_fraction=("absolute_gap", lambda x: np.mean(x < .01))).reset_index().to_csv(output / "DESCRIPTIVES.tsv", sep="\t", index=False)
    models.to_csv(output / "GENUS_MODELS.tsv", sep="\t", index=False)
    model_diagnostics.to_csv(output / "MODEL_DIAGNOSTICS.tsv", sep="\t", index=False)
    cross_validation.to_csv(output / "PASSPORT_CROSS_VALIDATION.tsv", sep="\t", index=False)
    pairs.to_csv(output / "MATCHED_HURWITZ_PAIRS.tsv", sep="\t", index=False)
    pair_summary.to_csv(output / "MATCHED_HURWITZ_SUMMARY.tsv", sep="\t", index=False)
    relations.to_csv(output / "BRAID_RELATION_DIAGNOSTICS.tsv", sep="\t", index=False)
    bounds.to_csv(output / "BIPARTITE_BOUNDS.tsv", sep="\t", index=False)
    traces.to_csv(output / "TRACE_MOMENT_ANALYSIS.tsv", sep="\t", index=False)
    convergence.to_csv(output / "LANCZOS_CONVERGENCE_ANALYSIS.tsv", sep="\t", index=False)
    high_reps.to_csv(output / "HIGH_REPLICATE_BASELINES.tsv", sep="\t", index=False)
    high_reps_summary.to_csv(output / "HIGH_REPLICATE_SUMMARY.tsv", sep="\t", index=False)
    switches.to_csv(output / "SWITCH_CHAIN_SENSITIVITY.tsv", sep="\t", index=False)
    recurring.to_csv(output / "RECURRING_EIGENVALUE_FAMILIES.tsv", sep="\t", index=False)
    make_figure(output, cases, models, traces, convergence, recurring)

    cv = cross_validation.pivot(index="outcome", columns="model", values="rmse")
    relation_text = "pending enhanced census"
    if len(relations):
        actual = relations[relations["source"] == "case"]
        relation_text = f"actual relation-defect fraction {actual.iloc[0]['fraction_with_relation_defect']:.3g}" if len(actual) else "actual diagnostics unavailable"
    trace_text = "Deep trace-moment run pending." if traces.empty else f"Deep trace analysis contains {len(traces)} baseline-by-order comparisons with Holm correction."
    report = f"""# Structural spectral follow-up analysis

This analysis treats all genus tests as exploratory observational comparisons. It describes {len(cases):,} quotient components from {census}; sampled S8 cases are excluded from genus models, matching, and passport validation. Regression uses passport-clustered standard errors and adjusts for group-degree, log graph size, loop density, parallel-edge density, and one generator cycle-density statistic when available. The design is checked for full rank. The leave-one-passport-out calculation tests prediction on unseen passports rather than fit to the existing catalog.

## Predeclared questions

1. Does genus improve relaxation-gap or absolute-gap prediction after structural adjustment?
2. Do matched Hurwitz components of different genus retain a systematic gap slope?
3. Do Hurwitz graphs have excess short closed walks relative to degree-, Schreier-, and generator-cycle-matched nulls?
4. Are the negative extremal eigenvalues and Lanczos estimates numerically stable?
5. Does the numerical value lambda2 near 2 + 2 sqrt(2) concentrate in recurring graph fingerprints, and can small examples be certified exactly?
6. How far from bipartite can the graph be certified, using loop bounds and numerical spectral estimates versus an explicit cut?

## Current checks

- Relation validation: {relation_text}.
- Leave-one-passport-out RMSE table: {cv.to_dict() if len(cv) else {}}.
- {trace_text}
- High-replicate component-level empirical tail tests and regular-null switch-chain sensitivity are reported when their bounded runs have completed.
- OLS normality and homoscedasticity diagnostics are reported; exploratory standard errors are clustered by group-degree and passport. Repeated isomorphic graph families can still reduce effective independent information.
- Matched genus slopes use disjoint case pairs and have no bootstrap confidence interval because cases recur across passports and graph families.
- Families of genus coefficients and trace-moment sign tests use Holm correction. Effect estimates and confidence intervals remain primary; adjusted p-values are secondary.
- Trace-moment comparisons average within group-degree/passport first, then use passport-cluster differences for the sign test and descriptive bootstrap interval.
- `BIPARTITE_BOUNDS.tsv` separates exact bipartiteness, a rigorous loop lower bound on frustration, a floating-point small-graph spectral estimate, and local-search feasible-cut upper bounds where the deep run evaluated them. The local-search bounds are heuristic; Lanczos minimum estimates are not used as certified lower bounds.

## Limitations

The empirical braid-aware control matches one real Hurwitz component to another, so it preserves braid relations but is not an independently sampled probability model. Cycle-matched Schreier nulls preserve every generator cycle type while deliberately breaking the braid relations. Neither comparison identifies a causal effect of genus, because genus is determined by the passport. Sampled S8 results are explicitly non-exhaustive.
The `2+2sqrt(2)` screen uses floating-point eigenvalues and Weisfeiler-Lehman fingerprints, neither of which is an exact certificate. The separately prepared exact-certificate script reruns eligible small cases, canonically identifies labeled graphs, and isolates roots of integer characteristic polynomials.
"""
    (output / "ANALYSIS.md").write_text(report, encoding="utf-8")
    (output / "FIGURE_ALT_TEXT.txt").write_text("Six panels compare relaxation and absolute gaps, the frequency of near-bipartite components, adjusted genus coefficients with 95% passport-clustered confidence intervals, normalized trace-moment differences from matched nulls, Lanczos convergence errors, and counts of graph fingerprints with a numerical eigenvalue near 2+2sqrt(2). Missing deep-run panels are labeled pending rather than interpolated.\n", encoding="utf-8")
    (output / "PROVENANCE.json").write_text(json.dumps({"seed": SEED, "census": str(census), "deep": str(deep) if deep else None, "convergence": str(convergence_root) if convergence_root else None, "high_reps": str(high_reps_root) if high_reps_root else None, "switch_runs": args.switch_run, "analysis": str(Path(__file__).resolve()), "multiple_testing": "Holm", "regression_covariance": "group-degree/passport clustered", "trace_uncertainty": "5000-replicate component bootstrap (descriptive)"}, indent=2) + "\n", encoding="utf-8")
    print(output / "ANALYSIS.md")


if __name__ == "__main__":
    main()
