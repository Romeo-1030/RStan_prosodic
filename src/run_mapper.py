import argparse
import numpy as np
import pandas as pd
import json
import os

from gtda.mapper import CubicalCover, make_mapper_pipeline
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS
from sklearn.base import BaseEstimator, TransformerMixin

# === Custom MDS Wrapper ===
class WrappedMDS(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=5, normalized_stress='auto', random_state=None):
        self.n_components = n_components
        self.normalized_stress = normalized_stress
        self.random_state = random_state
        self.mds = MDS(
            n_components=n_components, 
            dissimilarity='precomputed',
            normalized_stress=normalized_stress,
            random_state=random_state
        )

    def fit(self, X, y=None):
        self.mds.fit(X)
        return self

    def transform(self, X):
        return self.mds.transform(X)

    def fit_transform(self, X, y=None):
        return self.mds.fit_transform(X)

# === Build the Mapper pipeline ===
def build_pipeline(distance_matrix, n_intervals, overlap, eps, seed):
    filter_func = WrappedMDS(n_components=5, random_state=seed)
    cover = CubicalCover(n_intervals=n_intervals, overlap_frac=overlap)
    clusterer = DBSCAN(metric="precomputed", eps=eps).fit(distance_matrix)

    pipe = make_mapper_pipeline(
        filter_func=filter_func,
        cover=cover,
        clusterer=clusterer,
        verbose=False,
        n_jobs=1,
    )
    return pipe

# === Post-process the Mapper graph ===
def process_mapper_graph(mapper_graph):
    adjacency_matrix = np.array(mapper_graph.get_adjacency()).tolist()
    node_elements = [list(map(int, ne)) for ne in mapper_graph.vs["node_elements"]]

    to_remove = set()
    for i in range(len(node_elements)):
        for j in range(len(node_elements)):
            if i != j and set(node_elements[i]) < set(node_elements[j]):
                to_remove.add(i)
                break

    to_remove = sorted(to_remove, reverse=True)
    for idx in to_remove:
        del node_elements[idx]
        del adjacency_matrix[idx]
        for row in adjacency_matrix:
            del row[idx]

    return {
        "adjacency": adjacency_matrix,
        "points_in_vertex": node_elements
    }

# === Print node contents with labels ===
def print_nodes(export_data, column_names):
    for node_id, indices in enumerate(export_data["points_in_vertex"]):
        words = [column_names[i] for i in indices]
        print(f"Node {node_id}:")
        print(f"  Words: {words}")

# === Save result to a specific directory ===
def save_result(export_data, seed, n_intervals, overlap, eps, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    filename = (f"mapper_result_seed{seed}_int{n_intervals}" +
                f"_ov{str(overlap).replace('.', '_')}_eps{str(eps).replace('.', '_')}.json")
    output_file = os.path.join(out_dir, filename)
    with open(output_file, "w") as f:
        json.dump(export_data, f)
    print(f"\nSaved mapper result to: {output_file}")

# === Run everything ===
def run(seed, n_intervals, overlap, eps, out_dir):
    np.random.seed(seed)
    distance_matrix = pd.read_csv("../data/distance_matrix.csv").values
    distance_matrix_raw = pd.read_csv("../data/distance_matrix.csv")
    distance_matrix_sqrt = np.sqrt(distance_matrix)
    column_names = distance_matrix_raw.columns.tolist()

    pipe = build_pipeline(distance_matrix_sqrt, n_intervals, overlap, eps, seed)
    mapper_graph = pipe.fit_transform(distance_matrix_sqrt)
    export_data = process_mapper_graph(mapper_graph)

    print_nodes(export_data, column_names)
    save_result(export_data, seed, n_intervals, overlap, eps, out_dir)

# === Entry Point ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Mapper with configurable parameters")
    parser.add_argument("--seed", type=int, default=688, help="Random seed")
    parser.add_argument("--n_intervals", type=int, default=2, help="Number of intervals in the cover")
    parser.add_argument("--overlap", type=float, default=0.25, help="Overlap fraction in the cover")
    parser.add_argument("--eps", type=float, default=0.25, help="Epsilon value for DBSCAN")
    parser.add_argument("--out_dir", type=str, default=".", help="Directory to save the result")
    args = parser.parse_args()

    run(
        seed=args.seed,
        n_intervals=args.n_intervals,
        overlap=args.overlap,
        eps=args.eps,
        out_dir=args.out_dir
    )


# Note for Haoran:
# How to use it?
# python src/run_mapper.py --seed 688 --n_intervals 4 --overlap 0.3 --out_dir data/cluster_result

