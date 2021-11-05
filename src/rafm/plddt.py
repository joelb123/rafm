# -*- coding: utf-8 -*-
"""Find hydrated waters in structure."""
# standard library imports
from pathlib import Path
from typing import List
from typing import Optional
from typing import Tuple


# 3rd-party imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from statsdict import Stat

# module imports
from . import NAME
from . import VERSION
from .common import APP
from .common import STATS

# global constants
ATOM_REC: str = "ATOM"
ATOM_START_POS: int = 13
ATOM_STOP_POS: int = 15
ATOMS: Tuple[str, ...] = ("CA",)
B_FACTOR_START: int = 61
B_FACTOR_STOP: int = 65
REC_TYPE_START: int = 0
REC_TYPE_STOP: int = 4
DEFAULT_MIN_LENGTH = 20
DEFAULT_MIN_COUNT = 20
DEFAULT_PLDDT_LOWER_BOUND = 80
DEFAULT_PLDDT_UPPER_BOUND = 100
DEFAULT_PLDDT_CRITERION = 91.2
DEFAULT_LDDT_CRITERION = 0.8
DEFAULT_OUT_FILE_TYPE = "png"
DEFAULT_RESIDUE_CRITERION = 80
CRITERION_TYPE = "median"
MODULE_NAME = __name__.split(".")[0]
EMPTY_PATH = Path()


def bin_labels(bin_type, lower_bound, upper_bound=DEFAULT_PLDDT_UPPER_BOUND):
    """Create labels for bins of different quantities."""
    if upper_bound == DEFAULT_PLDDT_UPPER_BOUND:
        upper_label = ""
    else:
        upper_label = f"_{upper_bound}"
    return f"pLDDT{lower_bound}{upper_label}_{bin_type}"


def extract_b_factors(file_path: Path) -> List[float]:
    """Return an array of B factors from a PDB file specified by file_path."""
    if not file_path.exists():
        raise ValueError(f"PDB file {file_path} does not exist")
    with file_path.open("rU") as f:
        # parse b_factors out of PDB file
        b_factor_list = [
            float(rec[B_FACTOR_START:B_FACTOR_STOP])
            for rec in f.readlines()
            if (
                (len(rec) > B_FACTOR_STOP)
                and (rec[REC_TYPE_START:REC_TYPE_STOP] == ATOM_REC)
                and (rec[ATOM_START_POS:ATOM_STOP_POS] in ATOMS)
            )
        ]
    return b_factor_list


def compute_plddt_stats(
    file_path,
    lower_bound=DEFAULT_PLDDT_LOWER_BOUND,
    min_count=DEFAULT_MIN_COUNT,
    min_length=DEFAULT_MIN_LENGTH,
    upper_bound=DEFAULT_PLDDT_UPPER_BOUND,
):
    """Compute stats on pLDDTs for a PDB file specified by file_path."""
    plddts = np.array(extract_b_factors(file_path))
    n_pts = len(plddts)
    mean = np.NAN
    median = np.NAN
    n_trunc_obs = np.NAN
    trunc_mean = np.NAN
    trunc_median = np.NAN
    trunc_frac = np.NAN
    if n_pts >= min_length:
        mean = plddts.mean().round(2)
        median = np.median(plddts).round(2)
        obs = plddts[(plddts >= lower_bound) & (plddts <= upper_bound)]
        n_trunc_obs = len(obs)
        if len(obs) >= min_count:
            trunc_mean = obs.mean().round(2)
            trunc_median = np.median(obs).round(2)
            trunc_frac = round(n_trunc_obs / n_pts, 2)
    return (
        n_pts,
        mean,
        median,
        n_trunc_obs,
        trunc_frac,
        trunc_mean,
        trunc_median,
        str(file_path),
    )


@APP.command()
@STATS.auto_save_and_report
def plddt_stats(
    pdb_list: Optional[List[Path]],
    criterion: Optional[float] = DEFAULT_PLDDT_CRITERION,
    min_length: Optional[int] = DEFAULT_MIN_LENGTH,
    min_count: Optional[int] = DEFAULT_MIN_COUNT,
    lower_bound: Optional[int] = DEFAULT_PLDDT_LOWER_BOUND,
    upper_bound: Optional[int] = DEFAULT_PLDDT_UPPER_BOUND,
    file_stem: Optional[str] = MODULE_NAME,
) -> None:
    """Calculate stats on bounded pLDDTs from list of PDB model files."""
    results = []
    criterion_label = bin_labels(CRITERION_TYPE, lower_bound, upper_bound)
    stats_file_path = Path(f"{file_stem}_plddt_stats.tsv")
    n_models_in = len(pdb_list)
    STATS["models_in"] = Stat(n_models_in, desc="models read in")
    STATS["min_length"] = Stat(min_length, desc="minimum sequence length")
    STATS["min_count"] = Stat(
        min_length, desc="minimum # of selected residues"
    )
    STATS["plddt_lower_bound"] = Stat(
        lower_bound, desc="minimum bound per-residue"
    )
    STATS["plddt_upper_bound"] = Stat(
        upper_bound, desc="maximum bound per-residue"
    )
    STATS["plddt_criterion"] = Stat(
        criterion, desc=f"minimum bounded {CRITERION_TYPE} for selection"
    )
    for file_path in pdb_list:
        results.append(
            compute_plddt_stats(
                file_path,
                lower_bound=lower_bound,
                min_count=min_count,
                min_length=min_length,
                upper_bound=upper_bound,
            )
        )
    stats = pd.DataFrame(
        results,
        columns=(
            [
                "residues_in_pLDDT",
                "pLDDT_mean",
                "pLDDT_median",
                bin_labels("count", lower_bound, upper_bound),
                bin_labels("frac", lower_bound, upper_bound),
                bin_labels("mean", lower_bound, upper_bound),
                criterion_label,
                "file",
            ]
        ),
    )
    logger.info(f"Writing stats to {stats_file_path}")
    stats.sort_values(by=criterion_label, inplace=True, ascending=False)
    stats = stats.reset_index()
    stats.index.name = f"{NAME}-{VERSION}"
    del stats["index"]
    if (lower_bound == DEFAULT_PLDDT_LOWER_BOUND) and (
        upper_bound == DEFAULT_PLDDT_UPPER_BOUND
    ):
        file_col = stats["file"]
        del stats["file"]
        stats["LDDT_expect"] = (
            1.0
            - (
                (1.0 - (stats[criterion_label] / 100.0))
                * (1.0 - DEFAULT_LDDT_CRITERION)
                / (1.0 - DEFAULT_PLDDT_CRITERION / 100.0)
            )
        ).round(3)
        stats["passing"] = stats["LDDT_expect"] >= DEFAULT_LDDT_CRITERION
        stats["file"] = file_col
    stats.to_csv(stats_file_path, sep="\t")
    total_residues = int(stats["residues_in_pLDDT"].sum())
    STATS["total_residues"] = Stat(
        total_residues, desc="number of residues in all models"
    )
    selected_stats = stats[stats[criterion_label] >= criterion]
    n_models_selected = len(selected_stats)
    frac_models_selected = round(n_models_selected * 100.0 / n_models_in, 0)
    STATS["models_selected"] = Stat(
        n_models_selected,
        desc=f"models passing {criterion_label}>={criterion}",
    )
    STATS["model_selection_pct"] = Stat(
        frac_models_selected, desc="fraction of models passing, %"
    )
    selected_residues = int(selected_stats["residues_in_pLDDT"].sum())
    STATS["selected_residues"] = Stat(
        selected_residues, desc="residues in passing models"
    )


@APP.command()
def plddt_select_residues(
    criterion: Optional[float] = DEFAULT_PLDDT_CRITERION,
    min_length: Optional[int] = DEFAULT_MIN_LENGTH,
    min_count: Optional[int] = DEFAULT_MIN_COUNT,
    lower_bound: Optional[int] = DEFAULT_PLDDT_LOWER_BOUND,
    upper_bound: Optional[int] = DEFAULT_PLDDT_UPPER_BOUND,
    file_stem: Optional[str] = MODULE_NAME,
) -> None:
    """Select residues from files matching criterion."""
    stats_file_path = Path(f"{file_stem}_plddt_stats.tsv")
    stats = pd.read_csv(stats_file_path, sep="\t")
    criterion_label = bin_labels(CRITERION_TYPE, lower_bound, upper_bound)
    count_label = bin_labels("count", lower_bound, upper_bound)
    plddt_list = []
    residue_list = []
    file_list = []
    for unused_row_num, row in stats.iterrows():
        if (
            (row["residues_in_pLDDT"] >= min_length)
            and (row[criterion_label] >= criterion)
            and (row[count_label] >= min_count)
        ):
            plddts = extract_b_factors(Path(row["file"]))
            n_res = len(plddts)
            plddt_list += plddts
            residue_list += [i for i in range(n_res)]
            file_list += [row["file"]] * n_res
    df = pd.DataFrame(
        {"file": file_list, "residue": residue_list, "pLDDT": plddt_list}
    )
    out_file_path = Path(f"{file_stem}_plddt{lower_bound}_{criterion}.tsv")
    logger.info(f"Writing residue file {out_file_path}")
    df.to_csv(out_file_path, sep="\t")


@APP.command()
@STATS.auto_save_and_report
def plddt_plot_dists(
    criterion: Optional[float] = DEFAULT_PLDDT_CRITERION,
    lower_bound: Optional[int] = DEFAULT_PLDDT_LOWER_BOUND,
    upper_bound: Optional[int] = DEFAULT_PLDDT_UPPER_BOUND,
    file_stem: Optional[str] = MODULE_NAME,
    out_file_type: Optional[str] = DEFAULT_OUT_FILE_TYPE,
    residue_criterion: Optional[int] = DEFAULT_RESIDUE_CRITERION,
) -> None:
    """Plot histograms of per-model and per-residue pLDDT distributions."""
    stats_file_path = Path(f"{file_stem}_plddt_stats.tsv")
    res_file_path = Path(f"{file_stem}_plddt{lower_bound}_{criterion}.tsv")
    fig_file_path = Path(f"{file_stem}_dists.{out_file_type}")
    per_model = pd.read_csv(stats_file_path, sep="\t")
    per_model = per_model.fillna(0.0)
    criterion_label = bin_labels(CRITERION_TYPE, lower_bound, upper_bound)
    x_axis = r"$pLDDT$"
    plddt_col = bin_labels(CRITERION_TYPE, lower_bound, upper_bound)
    per_model[x_axis] = per_model[plddt_col]
    select_residues = pd.read_csv(res_file_path, sep="\t")
    n_select = len(per_model[per_model[criterion_label] >= criterion])
    n_models = len(per_model)
    n_select_residues = len(
        select_residues[select_residues["pLDDT"] >= residue_criterion]
    )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.ecdfplot(data=per_model, x=x_axis, ax=ax, color="darkblue")
    sns.ecdfplot(
        data=select_residues,
        x="pLDDT",
        ax=ax,
        linestyle="dashed",
        color="orange",
    )
    if upper_bound != DEFAULT_PLDDT_UPPER_BOUND:
        upper_bound_str = f"-{upper_bound}"
    else:
        upper_bound_str = ""
    ax.legend(
        labels=[
            rf'$pLDDT_{{{lower_bound}{upper_bound_str}}}$ of {n_models} "{file_stem}" models',
            rf"$pLDDT$ by residue of {n_select} passing models",
        ]
    )
    per_model_val = 100 - round(int(n_select * 100.0 / n_models), 0)
    h_offset = -1
    v_offset = -10
    ax.vlines(criterion, 0.0, per_model_val / 100.0, color="darkblue")
    ax.text(
        criterion - h_offset,
        (per_model_val - v_offset) / 200.0,
        (
            rf"$pLDDT_{{{lower_bound}{upper_bound_str}}}$"
            + f"\n = {criterion},\n"
            + f"{100-per_model_val}% pass\n"
            + "by model"
        ),
        color="darkblue",
    )
    per_residue_val = 100 - int(
        round(n_select_residues * 100.0 / len(select_residues), 0)
    )
    STATS["usable_residues_pct"] = Stat(
        100 - per_residue_val, desc=f"residues with LDDT > {residue_criterion}"
    )
    v_offset = 10
    ax.text(
        residue_criterion - h_offset,
        (per_residue_val - v_offset) / 200.0,
        (
            r"$pLDDT$"
            + f"\n= {residue_criterion},\n"
            + f"{100-per_residue_val}% pass\n"
            + "by residue"
        ),
        color="orange",
    )
    logger.info(f"Saving {fig_file_path}")
    sns.despine()
    plt.savefig(fig_file_path, dpi=300)
