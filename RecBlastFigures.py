
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os


def melt(df):
    """Melting the dataframe, returning a melted df, a list of species and a list of genes."""
    species_columns = [x for x in df.columns if x != 'gene_name']  # getting the columns with species names
    melted_df = pd.melt(df, id_vars=['gene_name'], value_vars=species_columns, var_name='Species',
                        value_name='Orthologues')  # melting
    melted_df.columns = ['Gene Name', 'Species', 'Orthologues']  # Changing column names
    # species list
    species = sorted(species_columns)
    # genes list
    genes = sorted(melted_df['Gene Name'].unique().tolist())
    return melted_df, species, genes


# receives melted_df
def create_swarmplot(df, path, title, colormap, genes, species):
    """
    The function creates a swarmplot using seaborn.
    :param df: pandas.DataFrame object
    :param path: The CSV file path.
    :param title: Title for the plot.
    :param colormap: Colormap
    :param genes: Ordered list of genes.
    :param species: Ordered list of species.
    :return:
    """
    print("Creating swarmplot for {}".format(path))
    output_path = os.path.dirname(path)
    output = join_folder(output_path, "%s_swarmplot.png" % title)
    fig = plt.figure(figsize=(16, 10), dpi=180)  # new figure
    sns.swarmplot(x='Gene Name', y='Orthologues', hue='Species', order=genes, hue_order=species, data=df,
                  palette=colormap)  # draw swarmplot
    plt.ylabel("# Orthologues")
    plt.xlabel("Gene Name")
    plt.ylim(0, )
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.savefig(output)  # saving figure as output
    plt.close()
    return output


# receives melted_df
def create_barplot(df, path, title, colormap, genes, species):
    """
    The function creates a bar plot using seaborn.
    :param df: pandas.DataFrame object
    :param path: The CSV file path.
    :param title: Title for the plot.
    :param colormap: Colormap
    :param genes: Ordered list of genes.
    :param species: Ordered list of species.
    :return:
    """
    print("Creating bar plot for {}".format(path))
    output_path = os.path.dirname(path)
    output = join_folder(output_path, "%s_barplot.png" % title)
    fig = plt.figure(figsize=(16, 10), dpi=180)
    sns.barplot(x='Gene Name', y='Orthologues', hue='Species', order=genes, hue_order=species, data=df,
                palette=colormap)
    plt.ylabel("#Orthologues")
    plt.xlabel("Gene Name")
    plt.ylim(0, )
    # plt.suptitle(title, fontsize=16)
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.savefig(output)
    plt.close()
    return output


# receives melted_df
def create_barplot_orthologues_by_species(df, path, title, colormap, genes, species):
    """
    The function creates a bar plot using seaborn.
    :param df: pandas.DataFrame object
    :param path: The CSV file path.
    :param title: Title for the plot.
    :param colormap: Colormap
    :param genes: Ordered list of genes.
    :param species: Ordered list of species.
    :return:
    """
    print("Creating barplot by species for {}".format(path))
    output_path = os.path.dirname(path)
    output = join_folder(output_path, "%s_barplot_byspecies.png" % title)
    fig = plt.figure(figsize=(16, 10), dpi=180)
    sns.barplot(x='Species', y='Orthologues', hue='Gene Name', data=df, order=species, hue_order=genes,
                palette=colormap)
    plt.ylabel("#Orthologues")
    plt.xlabel("Species")
    plt.ylim(0, )
    # plt.suptitle(title, fontsize=16)
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.savefig(output)
    plt.close()
    return output


def create_barplot_sum(df, path, title, colormap, species):
    """
    The function creates a bar plot of the sum of orthologues found using seaborn.
    :param df: pandas.DataFrame object
    :param path: The CSV file path.
    :param title: Title for the plot.
    :param colormap: Colormap
    :param species: Ordered list of species.
    :return:
    """
    print("Creating barplot of sum for {}".format(path))
    output_path = os.path.dirname(path)
    output = join_folder(output_path, "%s_barplot_sum.png" % title)
    fig = plt.figure(figsize=(16, 10), dpi=180)
    sns.barplot(x='Species', y='Orthologues', estimator=sum, ci=None, data=df, order=species, palette=colormap)
    plt.ylabel("#Orthologues")
    plt.xlabel("Species")
    plt.ylim(0, )
    # plt.suptitle(title, fontsize=16)
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.savefig(output)
    plt.close()
    return output


def create_heatmap(df, path, title, colormap):
    """
    Creates a heatmap for the Gene/Species orthologue data.
    :param title: The title of the figure
    :param df: a padnas.DataFrame object (not melted)
    :param path: Path for the input and output file
    :param colormap: Colormap
    :return:
    """
    print("Creating heatmap for {}".format(path))
    output_path = os.path.dirname(path)
    fig = plt.figure(figsize=(16, 10), dpi=180)
    plt.title(title, fontsize=16)
    sns.heatmap(df, annot=True, fmt="d", cmap=colormap)
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    output = join_folder(output_path, "%s_heatmap.png" % title)
    plt.savefig(output)
    plt.close()
    return output


def create_clustermap(df, path, title, colormap, col_cluster, dont_cluster):
    """
    Creates a cluster map for the Gene/Species orthologue data.
    :param dont_cluster: if True, skip clustering and return a blank image.
    :param col_cluster: if True, cluster the columns.
    :param title: The title of the figure
    :param df: a padnas.DataFrame object (not melted)
    :param path: Path for the input and output file
    :param colormap: Colormap
    :return:
    """
    print("Creating clustermap for {}".format(path))
    output_path = os.path.dirname(path)
    output = join_folder(output_path, "%s_clustermap.png" % title)
    fig = plt.figure(figsize=(16, 10), dpi=180)
    if not dont_cluster:  # if we want to cluster the columns (in case we have more than 2 columns)
        sns.clustermap(df, annot=True, col_cluster=col_cluster, fmt="d", cmap=colormap, linewidths=.5)
        # plt.suptitle(title, fontsize=16)
        plt.yticks(fontsize=10)
        plt.xticks(fontsize=10)
    plt.savefig(output)  # save the figure
    plt.close()
    return output


# generate heatmap and clustermap
def generate_visual_graphs(csv_rbh_output_filename, csv_strict_output_filename, csv_ns_output_filename):
    """
    The function generates heatmap + clustermap for the output data.
    :param csv_rbh_output_filename:
    :param csv_strict_output_filename:
    :param csv_ns_output_filename:
    :return:
    """
    # reading as data_frame (for heat/clustermaps)
    nonstrict_data = pd.read_csv(csv_ns_output_filename, index_col=0)
    strict_data = pd.read_csv(csv_strict_output_filename, index_col=0)
    rbh_data = pd.read_csv(csv_rbh_output_filename, index_col=0)

    # transpose data
    df_nonstrict = pd.DataFrame.transpose(nonstrict_data)
    df_strict = pd.DataFrame.transpose(strict_data)
    df_rbh = pd.DataFrame.transpose(rbh_data)

    # reading for other plots and melting them:
    melt_df_nonstrict, species_list, genes_list = melt(pd.read_csv(csv_ns_output_filename))
    melt_df_strict, species_list, genes_list = melt(pd.read_csv(csv_strict_output_filename))
    melt_df_rbh, species_list, genes_list = melt(pd.read_csv(csv_rbh_output_filename))
    print "Species list is: {}".format(species_list)
    print "Genes list is: {}".format(genes_list)

    # clustering enabler (( one is enough because all files contains the same amount of genes ))
    dont_cluster = False
    col_cluster = False
    if len(df_nonstrict.columns) > 2:
        col_cluster = True
    elif len(df_nonstrict.columns) < 2 or len(df_nonstrict) < 2:
        dont_cluster = True

    # visual outputs:
    viz_dict = dict()
    print("Creating heatmaps and clustermpaps")
    # non-strict
    viz_dict['non_strict_heatmap'] = create_heatmap(df_nonstrict, csv_ns_output_filename, 'non_strict', "BuGn")
    viz_dict['non_strict_clustermap'] = create_clustermap(df_nonstrict, csv_ns_output_filename, 'non_strict', "PuBu",
                                                          col_cluster, dont_cluster)
    viz_dict['non_strict_barplot'] = create_barplot(melt_df_nonstrict, csv_ns_output_filename, 'non_strict', "BuGn",
                                                    genes_list, species_list)
    viz_dict['non_strict_barplot_2'] = create_barplot_orthologues_by_species(melt_df_nonstrict, csv_ns_output_filename,
                                                                             'non_strict', "BuGn",
                                                                             genes_list, species_list)
    viz_dict['non_strict_swarmplot'] = create_swarmplot(melt_df_nonstrict, csv_ns_output_filename, 'non_strict', "BuGn",
                                                        genes_list, species_list)
    viz_dict['non_strict_barplotsum'] = create_barplot_sum(melt_df_nonstrict, csv_ns_output_filename,
                                                           'non_strict', "BuGn", species_list)
    # strict
    viz_dict['strict_heatmap'] = create_heatmap(df_strict, csv_strict_output_filename, 'strict', "Oranges")
    viz_dict['strict_clustermap'] = create_clustermap(df_strict, csv_strict_output_filename, 'strict', "YlOrRd",
                                                      col_cluster, dont_cluster)
    viz_dict['strict_barplot'] = create_barplot(melt_df_strict, csv_strict_output_filename, 'strict', "YlOrRd",
                                                genes_list, species_list)
    viz_dict['strict_barplot_2'] = create_barplot_orthologues_by_species(melt_df_strict, csv_strict_output_filename,
                                                                         'strict', "YlOrRd", genes_list, species_list)
    viz_dict['strict_swarmplot'] = create_swarmplot(melt_df_strict, csv_strict_output_filename,
                                                    'strict', "YlOrRd", genes_list, species_list)
    viz_dict['strict_barplotsum'] = create_barplot_sum(melt_df_strict, csv_strict_output_filename,
                                                       'strict', "YlOrRd", species_list)
    # RBH
    viz_dict['RBH_heatmap'] = create_heatmap(df_rbh, csv_rbh_output_filename, 'RBH', "YlGnBu")
    viz_dict['RBH_clustermap'] = create_clustermap(df_rbh, csv_rbh_output_filename, 'RBH', "YlGnBu", col_cluster,
                                                   dont_cluster)
    viz_dict['RBH_barplot'] = create_barplot(melt_df_rbh, csv_rbh_output_filename, 'RBH', "YlGnBu", genes_list,
                                             species_list)
    viz_dict['RBH_barplot_2'] = create_barplot_orthologues_by_species(melt_df_rbh, csv_rbh_output_filename, 'RBH',
                                                                      "YlGnBu", genes_list, species_list)
    viz_dict['RBH_swarmplot'] = create_swarmplot(melt_df_rbh, csv_rbh_output_filename, 'RBH', "YlGnBu",
                                                 genes_list, species_list)
    viz_dict['RBH_barplotsum'] = create_barplot_sum(melt_df_rbh, csv_rbh_output_filename, 'RBH', "YlGnBu",
                                                    species_list)

    return viz_dict


join_folder = os.path.join