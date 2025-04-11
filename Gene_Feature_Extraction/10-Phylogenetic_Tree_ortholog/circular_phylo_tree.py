#!/usr/bin/env python3
import ete3
print(ete3.__version__)


from ete3 import Tree, TreeStyle, TextFace

def circular_tree_visualization(newick_file, output_png="circular_tree.png"):
    """
    Reads a Newick format tree file, generates a circular (radial) tree visualization,
    and saves the output as a PNG image.

    Parameters:
        newick_file (str): Path to the Newick tree file (e.g., IQ-TREE .treefile).
        output_png (str): File name (or full path) where the PNG image should be saved.
    """
    try:
        # Load the tree from the provided file
        tree = Tree(newick_file, format=1)
    except Exception as e:
        print(f"[ERROR] Could not load tree from {newick_file}: {e}")
        return

    # Create a TreeStyle object for customization
    ts = TreeStyle()
    ts.mode = "c"               # Set mode to circular (radial)
    ts.show_leaf_name = True    # Show the names of the taxa at the leaves
    ts.title.add_face(TextFace("Circular Phylogenetic Tree", fsize=14), column=0)

    try:
        # Render the tree and save to a PNG file
        tree.render(output_png, w=800, units="px", tree_style=ts)
        print(f"[INFO] Circular tree visualization saved to: {output_png}")
    except Exception as e:
        print(f"[ERROR] Could not render the tree: {e}")

if __name__ == "__main__":
    # Specify the path to your Newick tree file. For example:
    tree_file = "phylo_tree_noBB.treefile"  # Replace with your actual path if needed.
    output_file = "circular_tree.png"         # The output image name or full path.
    
    circular_tree_visualization(tree_file, output_file)
