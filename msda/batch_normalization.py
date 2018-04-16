import pandas as pd
import numpy as np


def normalize_within_batch(df, samples, control='Bridge'):
    """Function normalizes sample in each batch by Bridge

    Parameters
    ----------
    df : pandas dataframe
        dataframe corresponding to mass spec results of one batch
    samples : list[str]
        list of sample names

    Returns
    -------
    df_norm_intra : pandas dataframe
       dataframe corresponding to mass spec results of one batch
       normalized by its bridge sample
    """

    df_norm_intra = df.copy()
    bridge_mean = float(df.loc[:, control].sum())
    true_samples = samples[:]
    true_samples.remove(control)

    for sample in true_samples:
        sample_mean = float(df.loc[:, sample].sum())
        nrm = bridge_mean / sample_mean
        # print(nrm)
        df_norm_intra.loc[:, sample] = nrm * df_norm_intra.loc[:, sample]
    return df_norm_intra


def normalize_between_batches(df, df_ref, samples, control='Bridge'):
    """Function normalizes samples across batches by ratio of Bridge samples

    Parameters
    ----------
    df : pandas dataframe
        dataframe corresponding to mass spec results of one batch
        that is already been normalzied wrt its bridge
    df_ref : pandas dataframe
        dataframe whose Bridge serves as primary normalizer
    samples : list[str]
        list of sample names

    Returns
    -------
    df_norm_inter : pandas dataframe
       dataframe corresponding to mass spec results of one batch
       normalized by its ratio of bridge samples across batches
    """

    true_samples = samples[:]
    true_samples.remove(control)
    refset_bridge_mean = float(df_ref.loc[:, 'Bridge'].sum())
    set_bridge_mean = float(df.loc[:, 'Bridge'].sum())
    nrm = refset_bridge_mean / set_bridge_mean
    df_norm_inter = df.copy()
    df_norm_inter.loc[:, true_samples] = nrm * df_norm_inter.loc[:,
                                                                 true_samples]
    return df_norm_inter


def normalize_mix(df, df_ref, samples, control='Bridge'):

    df_norm_mix = df.copy()
    ref_control = [s for s in df_ref.columns.tolist() if control in s][0]
    # true_samples = samples[:]
    # true_samples.remove(control)
    refset_bridge_mean = float(df_ref.loc[:, ref_control].sum())
    for sample in samples:
        sample_mean = float(df_norm_mix.loc[:, sample].sum())
        nrm = refset_bridge_mean / sample_mean
        df_norm_mix.loc[:, sample] = nrm * df_norm_mix.loc[:, sample]
    return df_norm_mix


def normalize_per_protein(df, df_refs, samples, control='Bridge'):
    df_pr = df.copy()
    # true_samples = samples[:]
    # true_samples.remove(control)
    proteins = df_pr.index.tolist()
    # missing_proteins = []

    for protein in proteins:
        ref_list = get_ref_list(df_refs, protein)
        df_ref = ref_list[0]
        ref_control = [s for s in df_ref.columns.tolist() if control in s][0]
        protein_in_refset = df_ref.loc[protein][ref_control]
        # try:
        #     protein_in_refset = df_ref[control][df_ref['Pro_Id']
        #                                         == protein].values[0]
        # except IndexError:
        #     missing_proteins.append(protein)
        pr_control = [s for s in df_pr.columns.tolist() if control in s][0]
        protein_in_set = df_pr[pr_control][df_pr.index == protein].values[0]
        if protein_in_set == 0:
            nrm = 1
        else:
            nrm = protein_in_refset / protein_in_set
        # index = df_pr[df_pr.Uniprot_Id == protein].index.values[0]
        #  print(protein, protein_in_set, protein_in_refset)
        df_pr.loc[protein, samples] = nrm * df_pr.loc[protein, samples]
    return df_pr


def get_ref_list(df_refs, protein):
    ref_list = []
    for dfr in df_refs:
        if protein in dfr.index.tolist():
            ref_list.append(dfr)
    return ref_list

