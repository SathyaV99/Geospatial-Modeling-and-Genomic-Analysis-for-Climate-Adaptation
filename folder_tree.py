import os

def generate_tree(path, prefix=''):
    tree_lines = []
    contents = sorted(os.listdir(path))
    files = [f for f in contents if os.path.isfile(os.path.join(path, f))]
    folders = [f for f in contents if os.path.isdir(os.path.join(path, f))]

    entries = folders + files
    for i, name in enumerate(entries):
        connector = "└── " if i == len(entries) - 1 else "├── "
        full_path = os.path.join(path, name)
        tree_lines.append(f"{prefix}{connector}{name}/" if os.path.isdir(full_path) else f"{prefix}{connector}{name}")

        if os.path.isdir(full_path):
            extension = "    " if i == len(entries) - 1 else "│   "
            tree_lines.extend(generate_tree(full_path, prefix + extension))

    return tree_lines

def save_tree_to_file(root_path, output_file):
    root_name = os.path.basename(os.path.normpath(root_path)) + '/'
    tree_structure = [root_name] + generate_tree(root_path)
    tree_text = "\n".join(tree_structure)

    # Print to console
    print(tree_text)

    # Save to .txt
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(tree_text)
    print(f"\n✅ Folder tree saved to: {output_file}")

# === USAGE ===
root_folder = r"D:\Documents\Python Stuff - Programming\AMOD Big Data research project\Genomic_and_Geographic_analysis_of_high_altitude_bovids"  # Change this to your target folder  # Change to your folder
output_txt = "root_folder.txt"  # Change this to your target folder  # Save output here

save_tree_to_file(root_folder, output_txt)
