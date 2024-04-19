from typing import List, Union

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from testsuite.algorithm_wrapper import AlgorithmWrapper


def seed_list_to_mask(seed_list: List[int], n: int) -> np.ndarray:
    """Convert list of seed nodes to a mask vector.

    Args:
    ----
    seed_list: List of nodes that are seeds.
    n: Number of entries in the returned array.

    Returns:
    -------
    Indicator array with 1s at seed indices.

    """
    train_seed_mask = np.zeros(n, dtype=int)
    train_seed_mask[seed_list] = 1
    return train_seed_mask


def qrw_score(
    G: nx.Graph,
    seed_list: List,
    t: float,
    H: csr_matrix,
    diag: Union[float, None] = None,
) -> np.ndarray:
    """Calculate quantum walk scores.

    Args:
    ----
    G: Graph that the walk occures on.
    seed_list: List of seed nodes.
    t: Time for which the walk lasts.
    H: Matrix to use as Hamiltonian.
    diag: How to set diagonals of the Hamiltonian.

    Returns:
    -------
    Array containing scores for each node in G.

    """
    n = G.number_of_nodes()
    n_seeds = len(seed_list)
    if isinstance(diag, (float, int)):
        # Construct sparse nxn matrix, with values `diag` at
        # entry (seed, seed), for each seed in seed_list:
        D = csr_matrix(([diag] * n_seeds, (seed_list, seed_list)), shape=(n, n))
        H += D
    # Trick for doing "quantum expm_multiply":
    Z = np.zeros((n, n_seeds), dtype=int)
    Z[seed_list, np.arange(n_seeds)] = 1  # each column corresponds to a seed
    res = expm_multiply(-1j * t * H, Z)
    res = np.abs(res) ** 2
    return res.sum(axis=1)


class CustomWrapper(AlgorithmWrapper):
    def run_algorithm(
        self,
        ggi_network,
        expression_data,
        phenotypes,
        seed_genes,
        p_values,
        indicator_matrix,
        prefix,
    ):
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
        p_values : dict of str: float
            P-values for all genes (keys are gene IDs).
        indicator_matrix : pd.DataFrame
            Indicator matrix obtained from expression data (indices are sample IDs, column names are gene IDs).
        prefix : str
            Prefix to be used for temporary files and directories.

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        mean_degree : float
            Mean degree of the result genes.
        """

        topn = 200
        result_genes = []
        G = ggi_network
        nl = np.arange(G.number_of_nodes())
        A = nx.adjacency_matrix(G, nodelist=nl)
        node_to_gid = {}
        gid_to_node = {}
        for node in G.nodes():
            node_to_gid[node] = G.nodes[node]["GeneID"]
            gid_to_node[G.nodes[node]["GeneID"]] = node
        seed_list = []
        missing = []
        for gid in seed_genes:
            try:
                seed_list.append(int(gid_to_node[gid]))
            except:
                missing.append(gid)
        print(f"{len(missing)} / {len(seed_genes)} missing seeds")
        scores = qrw_score(ggi_network, seed_list, H=A, t=0.45, diag=5)
        train_seed_mask = seed_list_to_mask(seed_list, G.number_of_nodes())
        test_mask = (1 - train_seed_mask).astype(bool)
        scores_test = scores[test_mask]
        ind = np.argpartition(scores_test, -topn)[-topn:]
        top_N = ind[np.argsort(scores_test[ind])]
        for i in range(topn):
            result_genes.append(node_to_gid[top_N[i]])  # reverse order?

        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
