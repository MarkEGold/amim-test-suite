from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd
from testsuite.utils import load_ggi_network
import numpy as np
import os
os.chdir("algorithms/")

class ClustEx2Wrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes):
        """Runs the algorithm.

        Parameters
        ----------
        ggi_network : nx.Graph
            Possibly randomized GGI network.
        expression_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        phenotypes : np.array, shape (n_samples,)
            Phenotype data (indices are sample IDs).

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        mean_degree : float
            Mean degree of the result genes.
        """

        # Write GGI network in format required by cosine.
        path_to_network = '../temp/cosine_ggi.txt'
        gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
        inv_map = {str(v): k for k, v in gene_ids.items()}
        expression_data = expression_data[list(inv_map.keys())]
        expression_data = expression_data.rename(columns=inv_map)
        cols = list(expression_data.columns)
        cols = sorted(cols)
        expression_data = expression_data[cols]
        with open(path_to_network, 'w') as edge_list_file:
            for u, v in ggi_network.edges():
                edge_list_file.write(f'{u}{","}{v}\n')

        expression_data.loc[expression_data.index[phenotypes == 1]].to_csv('../temp/cosine_expr1.txt')
        expression_data.loc[expression_data.index[phenotypes == 0]].to_csv('../temp/cosine_expr2.txt')

        cosine = 'cd ../algorithms/cosine/; Rscript cosine.R'
        subprocess.call(cosine, shell=True, stdout=subprocess.PIPE)

        # Read the results.
        result_genes = []
        path_to_output = '../temp/cosine_output.txt'
        with open(path_to_output, 'r') as results:
            for line in results:
                result_genes.append(gene_ids[int(line.strip())])

        # Delete temporary data.
        subprocess.call('rm ../temp/cosine_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)