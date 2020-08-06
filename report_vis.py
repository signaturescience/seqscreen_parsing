#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 13:51:33 2020

@author: Student
"""
import argparse
import os
import pandas as pd
import altair as alt
from Bio import Entrez

pd.set_option("display.max_colwidth", 10000)

def get_strname(dframe):
    """
    Extracts one common string name for each group.

    Parameters
    ----------
    dframe : dataframe representation of the groups file.

    Returns
    -------
    a dataframe with the column "grlabel" added, containing standardized labels
    for each group.

    """
    dfs = {}
    final = {}
    master = {}

    for grp in list(dframe["Group"].unique()):
        dfs[grp] = list(dframe.loc[dframe["Group"] == grp, "Sample"].str.split("_"))
        final[grp] = []
        for lst in dfs[grp]:
            for word in lst:
                if word in master.keys():
                    master[word] += 1
                else:
                    master[word] = 1

    for grp, lst in dfs.items():
        for lst2 in lst:
            final[grp] = " ".join([w for w in lst2 if (master[w] >= len(
                dframe.loc[dframe["Group"] == grp].index) and w.lower() != "covid19")])
    dframe["grlabel"] = dframe["Group"].copy()
    return dframe.replace({"grlabel":final})

def get_grname(dframe):
    """
    Extracts a string name in common from groups doc.

    Parameters
    ----------
    dframe : dataframe representation of the groups file.

    Returns
    -------
    labels : a dataframe containing group number, file/sample name and group label

    """
    labels = {}
    for grp in list(dframe["group"].unique()):
        temp = dframe.loc[dframe["group"] == grp, "group":"publication"]
        labels[grp] = temp["publication"].unique()[0].split("_")[0]
        temp2 = temp["type"].str.split("_", expand=True)
        labels[grp] = " ".join([labels[grp], temp2.iloc[:, -1].unique()[0]])

    dframe.insert(1, "grlabel", dframe["group"].copy().replace({"grlabel":labels}))
    return pd.concat([dframe["Group"], dframe["Sample"], dframe["grlabel"]], 1)

def kraken_altair(dirname, groups, taxids):
    """
    Parameters
    ----------
    dirname : path to kraken report directory (string).
    groups : path to document specifying groups of sample files in dirname.
    taxids : list of taxids (integers) searched for in the input reports

    Returns
    -------
    fin : dataframe compatible with altair visualisation.

    """
    #ignores hidden files in input directory
    inlist = [f for f in os.listdir(dirname) if not f.startswith('.')]
    frames = []
    grp = get_strname(pd.read_csv(groups, sep="\t", index_col=False))
    nin = []

    for file in inlist:
        fname = file.split("_")[0] #label for each bar
        dfr = pd.read_csv(f"{dirname}/{file}", sep="\t", skiprows=0)
        dfr.columns = range(dfr.shape[1]) #resets column names

        if taxids is None: #breaks down "other" based on threshold
            temp0 = dfr.where(dfr[3].str.strip() == "S").dropna()
            temp = pd.concat([temp0[0], temp0[5].str.strip()], 1)
            temp.loc[len(temp)] = [100 - float(dfr.iloc[0, 0]), "unclassified"]

        else:
            frs = [pd.concat([dfr[0].where(dfr[4] == int(tid)).astype(float),
                              dfr[5].str.strip().where(dfr[4] == int(tid)).dropna()],
                             1).dropna() for tid in taxids]
            temp = pd.concat(frs)
            temp.loc[len(temp)] = [100 - temp[0].sum(), "other"]

        temp["file"] = fname
        ind = grp.index[grp["Sample"].str.contains(fname)]

        if not ind.empty:
            ind = ind.values[0]
            temp["group"] = grp.iat[ind, 0]
            temp["grlabel"] = grp.iat[ind, len(grp.columns) -1].split(f"_{fname}")[0]
            frames.append(temp)
        else:
            nin.append(fname)
    fin = pd.concat(frames, 0, ignore_index=True)
    fin = fin.rename(columns={0:"taxon abundances", 5:"species"})
    if nin != []:
        print(f"Warning: the following files can't be found in '{groups}': "
              + f"{(str(nin).split('[')[-1]).split(']')[0]}")
    return fin

def get_rank(taxid, rank, email):
    """

    Parameters
    ----------
    taxid : taxonomic identifier.
    rank : str representing desired rank (must be valid in the NCBI database).

    Returns
    -------
    string representing classification of the same organism at desired rank.

    """
    Entrez.email = email
    search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")
    record = Entrez.read(search)[0]
    if record["Rank"] == rank:
        return record['ScientificName']
    else:
        for dic in record["LineageEx"]:
            if dic["Rank"] == rank:
                return dic['ScientificName']

def seqscreen_altair(dirname, rank, email):
    """

    Parameters
    ----------
    dirname : path to seqscreen report directory (string).

    Returns
    -------
    dataframe compatible with altair visualisation.

    """
    inlist = [f for f in os.listdir(dirname) if not f.startswith('.')]
    dframes = []
    for file in inlist:
        path = f"{dirname}/{file}"
        print("path", path)
        temp = pd.read_csv(path, sep="\t")
        if rank is not None:
            tids = temp["taxid"].to_frame().apply(get_rank, 1, args=[rank, email])
        else:
            tids = temp["organism"].rename(0)
        tids = tids.value_counts(normalize=True).multiply(100)
        tids = tids.reset_index(
            ).rename(columns={"index":"species", 0:"taxon abundances"})
        tids["file"] = (file.split("_")[0]).split(".")[0]
        tids["grlabel"] = " "
        dframes.append(tids)
    fin = pd.concat(dframes, 0, ignore_index=True)
    print(fin)
    return fin

def stacked_bar(dframe, destname, kraken):
    """

    Parameters
    ----------
    dframe : dataframe compatible with altair.
    destname : path to where chart will be saved.
    kraken: whether the input data is a Kraken or Seqscreen report (True for
    Kraken reports, False for Seqscreen reports)

    Generates interactive stacked and grouped bar chart displaying taxon
    abundances across groups and samples.

    """
    sel1 = alt.selection_multi(fields=["species"], bind="legend")
    sel2 = alt.selection_multi(fields=["species"])
    scale = alt.Scale(domain=(0, 100))

    chart = alt.Chart(dframe).mark_bar().encode(
        x=alt.X("file" if kraken else "sum(taxon abundances):Q",
                **{"title":None} if kraken else {}),
        y=alt.Y("sum(taxon abundances):Q" if kraken else "file",
                **{"scale":scale} if kraken else {}),
        column=alt.Column('grlabel:N', title=None, spacing=50),
        order=alt.Order('species:N', sort='ascending'),
        color=alt.condition(
            (sel1 | sel2),
            alt.Color("species:N", scale=alt.Scale(scheme=("set1" if kraken else "sinebow"))),
            alt.value('lightgray'))
        ).resolve_scale(x='independent').add_selection(sel1, sel2).properties(
            width=(0 if kraken else 1000), height=(0 if kraken else 500)
            ).configure_legend(labelLimit=0,
                               symbolLimit=len(dframe.index)+1,
                               ).configure_axis(domain=False)
    chart.save(f'{destname}_chart.html')

def main():
    """
    makes bar charts for kraken and seqscreen outputs from command line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dirname", type=str,
        help="path to directory containing kraken2 or seqscreen reports")

    parser.add_argument("-k", "--kraken", type=str, nargs="+")
    parser.add_argument("-s", "--seqscreen", type=str, default=None, nargs="+")
    args = parser.parse_args()

    if args.seqscreen:
        sdf = seqscreen_altair(args.dirname, args.seqscreen[0], args.seqscreen[1])
        stacked_bar(sdf, args.dirname, False)

    elif args.kraken:
        groups = args.kraken[0]
        args.kraken.pop(0)
        kdf = kraken_altair(args.dirname, groups, args.kraken)
        stacked_bar(kdf, args.dirname, True)

main()
