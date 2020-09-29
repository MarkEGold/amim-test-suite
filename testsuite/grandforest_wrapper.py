from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd


class GrandForestWrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes, gene_scores, indicator_matrix):
        """Runs the algorithm.

        Parameters
        ----------
        ggi_network : nx.Graph
            Possibly randomized GGI network.
        expression_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        phenotypes : np.array, shape (n_samples,)
            Phenotype data (indices are sample IDs).
        seed_genes : list of str
            Seed genes (entries are gene IDs).
        gene_scores : dict of str: float
            Scores for all genes (keys are gene IDs).
        indicator_matrix : pd.DataFrame
            Indicator matrix obtained from expression data (indices are sample IDs, column names are gene IDs).

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        mean_degree : float
            Mean degree of the result genes.
        """

        # Write GGI network in format required by ClustEx2.
        path_to_network = '../temp/gf_ggi.txt'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_to_network, '\t', 'source\ttarget')

        gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
        ppi_genes = list(gene_ids.values())
        # ppi_genes = [int(x) for x in ppi_genes]
        # expression_data = expression_data[ppi_genes]
        expression_data["phenotype"] = phenotypes

        expression_data.to_csv('../temp/gf_expr.txt', index = False)

        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_to_network, '\t', 'source\ttarget')

        # Run GF
        grandforest = 'cd ../algorithms/grand_forest/; Rscript grandforest.R'
        subprocess.call(grandforest, shell = True, stdout=subprocess.PIPE)

        # Read the results.
        result_genes = []
        path_to_output = '../temp/gf_output.txt'
        with open(path_to_output, 'r') as results:
            for line in results:
                result_genes.append(line.strip())

        # Delete temporary data.
        subprocess.call('rm ../temp/gf_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
