import requests
import networkx as nx
import matplotlib.pyplot as plt
import csv
import os

G = None
data_global = []
main_proteins_global = []

def get_string_interactions(protein_names, species=9606, score_threshold=400):
    """
    Fetch protein-protein interactions from STRING database.
    
    Parameters:
        protein_names: list of str
        species: int, NCBI taxonomy ID
        score_threshold: int, minimum confidence score (0-1000)
    
    Returns:
        list of interactions
    """
    all_data = []
    for protein_name in protein_names:
        url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": protein_name,
            "species": species,
            "required_score": score_threshold
        }
        try:
            response = requests.get(url, params=params)
            data = response.json()
            if isinstance(data, list):
                all_data.extend(data)
        except Exception as e:
            print(f"Warning: failed to fetch {protein_name}: {e}")
            continue
    return all_data

def get_layout(G, layout="spring"):
    if layout == "spring":
        return nx.spring_layout(G, seed=42)
    elif layout == "circular":
        return nx.circular_layout(G)
    elif layout == "shell":
        return nx.shell_layout(G)
    else:
        return nx.spring_layout(G, seed=42)

def build_network(data, protein_names, show_plot=False, layout="spring"):
    """
    Build and optionally plot a protein interaction network.
    """
    global G, data_global, main_proteins_global
    G = nx.Graph()
    main_proteins = set([p.strip() for p in protein_names])
    main_proteins_global = main_proteins
    data_global = data

    for interaction in data:
        a = interaction.get('preferredName_A')
        b = interaction.get('preferredName_B')
        if a and b:
            G.add_node(a)
            G.add_node(b)
            G.add_edge(a, b, weight=interaction.get('score', 0))

    if show_plot:
        node_colors = ['red' if node in main_proteins else 'lightblue' for node in G.nodes()]
        pos = get_layout(G, layout)
        plt.figure(figsize=(10,8))
        nx.draw(G, pos, with_labels=True, node_color=node_colors,
                node_size=1500, font_size=10, edge_color='gray')
        plt.show()

    return G

def generate_report(G, protein_names, data, save_path="protein_network_report.txt"):
    """
    Generate a text report of the protein network.
    """
    main_proteins = set([p.strip() for p in protein_names])
    with open(save_path, "w") as f:
        f.write("=== Protein Network Report ===\n\n")
        f.write(f"Main proteins: {', '.join(main_proteins)}\n")
        f.write(f"Number of nodes: {G.number_of_nodes()}\n")
        f.write(f"Number of edges: {G.number_of_edges()}\n\n")

        degrees = dict(G.degree())
        if degrees:
            max_deg_node = max(degrees, key=degrees.get)
            f.write(f"Most connected protein: {max_deg_node} ({degrees[max_deg_node]} connections)\n\n")

        sorted_edges = sorted(data, key=lambda x: x.get('score', 0), reverse=True)[:5]
        f.write("Top 5 interactions:\n")
        for edge in sorted_edges:
            f.write(f"{edge.get('preferredName_A')} - {edge.get('preferredName_B')}, score: {edge.get('score')}\n")

        scores = [edge.get('score', 0) for edge in data]
        if scores:
            f.write(f"\nMax score: {max(scores)}\n")
            f.write(f"Min score: {min(scores)}\n")
            f.write(f"Average score: {sum(scores)/len(scores):.2f}\n")

    print(f"Report saved to {save_path}")
    os.startfile(save_path)

