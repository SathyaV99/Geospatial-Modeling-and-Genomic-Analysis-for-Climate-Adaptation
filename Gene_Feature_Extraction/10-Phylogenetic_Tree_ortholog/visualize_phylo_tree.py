import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo

def visualize_tree(tree_file, output_png="tree_visualization.png"):
    """
    Reads a phylogenetic tree in Newick format and visualizes it using Biopython.

    Parameters:
    -----------
    tree_file : str
        Path to the input tree file (Newick format, e.g., IQ-TREE .treefile).
    output_png : str
        File name (or path) for saving the plotted tree as a PNG image.
        If None, no figure will be saved to file.

    Notes:
    ------
    - This function displays a quick tree plot. For more advanced visualizations,
      consider specialized tools like FigTree or iTOL.
    - If running this code in a script (not an interactive environment),
      you may not see a pop-up window, so we rely on saving the PNG.
    """

    # Read the tree in Newick format
    tree = Phylo.read(tree_file, "newick")

    # Create a matplotlib figure
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # Draw the tree using Biopython
    Phylo.draw(
        tree,
        do_show=False,  # We'll manage the figure show manually
        axes=ax
    )

    # Optionally save to a PNG file
    if output_png:
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        print(f"[INFO] Tree visualization saved to {output_png}")
    else:
        print("[INFO] No output file specified. Displaying plot interactively.")
        plt.show()

    # Clean up the figure to avoid overlap if re-run in the same session
    plt.close(fig)


# Path to your .treefile from IQ-TREE
tree_file = "phylo_tree_noBB.treefile"

# Visualize and save as PNG
visualize_tree(tree_file, "my_bovid_tree.png")
