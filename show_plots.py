from os import listdir
from typing import List

import matplotlib.pyplot as plt
import mygene
import pandas as pd
import seaborn as sns


def flatten(genes: List) -> List:
    return [item for sublist in genes for item in sublist]


mypath = "results/"
all_results = listdir(mypath)


def get_pathways(condition: str) -> List[str]:
    if condition == "ALS":
        return ["hsa05014"]
    if condition == "LC":
        return ["hsa05223"]
    if condition == "UC":
        return ["hsa04060", "hsa04630", "hsa05321"]
    if condition == "HD":
        return ["hsa05016"]
    if condition == "CD":
        return ["hsa04621", "hsa04060", "hsa04630", "hsa05321", "hsa04140"]
    return []


new_files = [
    "HPRD_ORIGINAL_CUSTOM.csv",
    "HPRD_REWIRED_CUSTOM.csv",
    "HPRD_ORIGINAL_DIAMOND.csv",
    "HPRD_REWIRED_DIAMOND.csv",
    "HPRD_UNIFORM_DIAMOND.csv",
    "HPRD_UNIFORM_CUSTOM.csv",
    "HPRD_EXPECTED_DEGREE_CUSTOM.csv",
    "HPRD_EXPECTED_DEGREE_DIAMOND.csv",
]
new_files += [
    "CUSTOM_ORIGINAL_CUSTOM.csv",
    "CUSTOM_ORIGINAL_DIAMOND.csv",
    "CUSTOM_REWIRED_CUSTOM.csv",
    "CUSTOM_REWIRED_DIAMOND.csv",
    "CUSTOM_UNIFORM_CUSTOM.csv",
    "CUSTOM_UNIFORM_DIAMOND.csv",
    "CUSTOM_EXPECTED_DEGREE_CUSTOM.csv",
    "CUSTOM_EXPECTED_DEGREE_DIAMOND.csv",
]

print(new_files)
for i in range(len(new_files)):
    if new_files[i].endswith(".csv"):
        if i == 0:
            results = pd.read_csv(mypath + new_files[i])
        else:
            temp = pd.read_csv(mypath + new_files[i])
            results = pd.concat([results, temp], axis=0, ignore_index=True)

            COND = "GSE3790"
results = results[results.condition_name == COND].reset_index(drop=True)
results = results.replace({"algorithm_name": {"CUSTOM": "QA", "DIAMOND": "DIA"}})
results = results.replace(
    {
        "condition_name": {
            "GSE112680": "ALS",  # 198 seeds
            "GSE30219": "LC",  # 698
            "GSE75214": "UC",  # 1881
            "GSE75214_cd": "CD",  # 446
            "GSE3790": "HD",  # 634
        },
    },
)

all_genes = list(results.result_genes)
genes = []
for i in range(len(all_genes)):
    try:
        genes.append(all_genes[i].split(","))
    except:
        pass


genes = list(set(flatten(genes)))
mg = mygene.MyGeneInfo()
out = mg.querymany(
    genes,
    scopes="entrezgene",
    fields="symbol",
    species="human",
    verbose=False,
)
mapping = dict()
for line in out:
    try:
        mapping[line["query"]] = line["symbol"]
    except KeyError:
        print("{0} was not mapped to any gene name".format(line["query"]))
        mapping[line["query"]] = line["query"]


disgenet = pd.read_csv("data/networks/disgenet.csv")
dis_ids = {
    "ALS": "C0002736",
    "LC": "C1737250",
    "UC": "C0009324",
    "HD": "C0020179",
    "CD": "C0021390",
}
disgenet = disgenet[disgenet["disease_id"].isin(list(dis_ids.values()))]
overlaps = []
for i in results.index:
    try:
        gs = results.result_genes[i]

        condition = results["condition_name"][i]
        dis = disgenet[disgenet.disease_id == dis_ids[condition]]
        dis_genes = list(dis.gene)
        dis_genes = [str(x) for x in dis_genes]
        gs = gs.split(",")
        oc = len(set(gs).intersection(set(dis_genes))) / min(
            len(set(gs)),
            len(set(dis_genes)),
        )
        overlaps.append(oc)
    except:
        overlaps.append(0)

results["disgenet_overlap"] = overlaps

print(results)
sort_order = [
    "ORIGINAL",
    "REWIRED",
    "EXPECTED_DEGREE",
    "UNIFORM",
]
results["network_generator_name"] = pd.Categorical(
    results["network_generator_name"],
    categories=sort_order,
    ordered=True,
)
results = results.sort_values(by="network_generator_name")
results = results.replace(
    {
        "network_generator_name": {
            "ORIGINAL": "Original",
            "REWIRED": "Rewired",
            "UNIFORM": "Uniform",
            "EXPECTED_DEGREE": "Expected Degree",
        },
    },
)
results = results.rename(
    columns={
        "network_generator_name": "Generator",
        "mean_mutual_information": "Mean Mutual Information",
        "algorithm_name": "Model",
    },
)


def plot1() -> None:
    fig, axes = plt.subplots(1, 1)
    past = sns.color_palette()
    palette = {}
    palette["QA"] = "#54B6B8"
    palette["DIA"] = past[1]
    kv = {
        "x": "Generator",
        "data": results,
        # "hue": "condition_name",
        "hue": "Model",
        "palette": palette,
    }

    axes.set_title(f"Mean mutual information with phenotype for {COND}")
    sns.boxplot(ax=axes, y="Mean Mutual Information", **kv)
    plt.tight_layout()
    if COND:
        plt.savefig(f"img/basic/{COND}.png", dpi=600)
        plt.savefig(f"img/basic/{COND}.pdf", dpi=600)
    else:
        plt.savefig("img/basic/ALL.png", dpi=600)
        plt.savefig("img/basic/ALL.pdf", dpi=600)


def plot2() -> None:
    fig, axes = plt.subplots(2, 1)
    fig.subplots_adjust(hspace=0.8, wspace=0.4)
    past = sns.color_palette()
    palette = {}
    palette["QA"] = "#54B6B8"
    palette["DIA"] = past[1]
    kv = {
        "x": "Generator",
        "data": results,
        # "hue": "condition_name",
        "hue": "Model",
        "palette": palette,
    }

    sns.boxplot(ax=axes[0], y="Mean Mutual Information", **kv)
    axes[0].set_title(
        r"$\bf{" + "a)" + "}$" + " Mean mutual information with the phenotype",
        loc="left",
    )
    axes[0].legend_.remove()

    sns.boxplot(ax=axes[1], y="disgenet_overlap", **kv)
    axes[1].set_title(
        r"$\bf{" + "b)" + "}$" + " Overlap with DisGeNET disease genes",
        loc="left",
    )
    axes[1].legend_.remove()
    axes[1].set(xlabel="")

    plt.tight_layout()
    if COND:
        plt.savefig(f"img/basic/{COND}--2.png", dpi=600)
        plt.savefig(f"img/basic/{COND}--2.pdf", dpi=600)
    else:
        plt.savefig(f"img/basic/ALL--2.png", dpi=600)
        plt.savefig(f"img/basic/ALL--2.pdf", dpi=600)


plot1()
plot2()
