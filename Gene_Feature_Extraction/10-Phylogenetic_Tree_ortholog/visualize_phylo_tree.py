import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo

def visualize_tree(tree_file, output_png="tree_visualization.png"):


    # Read the tree in Newick format
    tree = Phylo.read(tree_file, "newick")

    # figure generate
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # draw phylograpgic tree
    Phylo.draw(
        tree,
        do_show=False,  
        axes=ax
    )

    if output_png:
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        print(f"[INFO] Tree visualization saved to {output_png}")
    else:
        print("[INFO] No output file specified. Displaying plot interactively.")
        plt.show()
    plt.close(fig)


tree_file = "phylo_tree_noBB.treefile"
visualize_tree(tree_file, "my_bovid_tree.png")
