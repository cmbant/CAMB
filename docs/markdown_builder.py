#!/usr/bin/env python
"""
This script builds Sphinx documentation in Markdown format and combines it into a single file
for use as context with Large Language Models (LLMs).

It can be used:
1. As a pre-build step in ReadTheDocs
2. Locally to generate markdown documentation
3. In CI/CD pipelines

Usage:
    python markdown_builder.py [--exclude file1,file2,...] [--output output_file] [--no-install]

Options:
    --exclude: Comma-separated list of files to exclude (without .md extension)
    --output: Output file path
"""

import os
import sys
import subprocess
import argparse
import glob
import traceback
import shutil


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Build Sphinx documentation in Markdown format for LLM context."
    )
    parser.add_argument(
        "--exclude",
        type=str,
        default="",
        help="Comma-separated list of files to exclude (without .md extension)",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file path",
    )
    return parser.parse_args()


def build_markdown_docs():
    """Build the documentation in Markdown format."""
    print("Building documentation in Markdown format...")
    build_dir = "docs/_build/markdown"

    # Create build directory if it doesn't exist
    os.makedirs(build_dir, exist_ok=True)

    # Create a temporary conf.py with intersphinx_disabled_reftypes setting
    temp_conf_dir = os.path.join(os.path.dirname(build_dir), "temp_conf")
    os.makedirs(temp_conf_dir, exist_ok=True)
    temp_conf_path = os.path.join(temp_conf_dir, "conf.py")
    with open("docs/source/conf.py", "r", encoding="utf-8") as f:
        conf_content = f.read()

    # Disable intersphinx extension for markdown build
    conf_content = conf_content.replace("'sphinx.ext.intersphinx', ", "")

    with open(temp_conf_path, "w", encoding="utf-8") as f:
        f.write(conf_content)

    try:
        # Run sphinx-build with the temporary conf.py
        result = subprocess.run(
            [
                "sphinx-build",
                "-b",
                "markdown",
                "-c",
                temp_conf_dir,
                "docs/source",
                build_dir,
            ],
            check=False,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"Warning: sphinx-build returned non-zero exit code: {result.returncode}")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            # Continue anyway as we might still have generated some markdown files
    finally:
        # Clean up the temporary conf directory
        if os.path.exists(temp_conf_dir):
            shutil.rmtree(temp_conf_dir)

    return build_dir


def extract_toctree_order(index_rst_path="docs/source/index.rst"):
    """Extract the order of files from the index.rst toctree."""
    try:
        with open(index_rst_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Find all toctree sections
        toctree_sections = []
        current_section = []
        in_toctree = False

        for line in content.split('\n'):
            line = line.strip()
            if '.. toctree::' in line:
                in_toctree = True
                current_section = []
            elif in_toctree:
                if line and not line.startswith(':'):
                    # This is a document reference in the toctree
                    current_section.append(line)
                elif not line and current_section:
                    # Empty line after entries - end of this toctree
                    toctree_sections.append(current_section)
                    in_toctree = False

        # Add the last section if we're still in a toctree at the end
        if in_toctree and current_section:
            toctree_sections.append(current_section)

        # Flatten the list of sections
        ordered_files = []
        for section in toctree_sections:
            ordered_files.extend(section)

        return ordered_files
    except Exception as e:
        print(f"Warning: Could not extract toctree order from index.rst: {e}")
        return []


def combine_markdown_files(build_dir, exclude_files, output_file):
    """Combine Markdown files into a single file with improved structure."""
    print(f"Combining Markdown files into {output_file}...")

    # Get all markdown files
    md_files = sorted(glob.glob(os.path.join(build_dir, "*.md")))

    if not md_files:
        print(f"Error: No markdown files found in {build_dir}")
        return False

    # Convert exclude_files to a set for faster lookup
    exclude_set = {f"{name.strip()}.md" for name in exclude_files if name.strip()}

    # Print excluded files for debugging
    if exclude_set:
        print(f"Excluding the following files: {', '.join(exclude_set)}")

    # Filter out excluded files
    filtered_files = [f for f in md_files if os.path.basename(f) not in exclude_set]

    # Check if any files were actually excluded
    excluded_count = len(md_files) - len(filtered_files)
    if excluded_count > 0:
        print(f"Successfully excluded {excluded_count} file(s)")
    else:
        print("Note: No files were excluded.")

    # Ensure index.md comes first if it exists
    index_file = None
    other_files = []

    for file_path in filtered_files:
        if os.path.basename(file_path).lower() == "index.md":
            index_file = file_path
        else:
            other_files.append(file_path)

    # Get the order from index.rst
    ordered_docs = extract_toctree_order()

    # Create a mapping of filenames to their paths
    filename_to_path = {os.path.splitext(os.path.basename(f))[0]: f for f in other_files}

    # Sort other_files based on the order in index.rst
    sorted_other_files = []

    # First add files in the order they appear in index.rst
    ordered_count = 0
    for doc in ordered_docs:
        if doc in filename_to_path:
            sorted_other_files.append(filename_to_path[doc])
            # Remove from the dictionary to mark as processed
            del filename_to_path[doc]
            ordered_count += 1

    # Then add any remaining files that weren't in the toctree
    remaining_files = list(filename_to_path.values())
    remaining_files.sort()  # Sort alphabetically
    sorted_other_files.extend(remaining_files)

    if ordered_count > 0:
        print(f"Sorted {ordered_count} files according to their order in index.rst")
    if remaining_files:
        print(f"Added {len(remaining_files)} additional files not found in index.rst (sorted alphabetically)")

    # Reorder files with index first, then the rest in the determined order
    if index_file:
        filtered_files = [index_file] + sorted_other_files
        print("Placing index.md first in the combined document")
    else:
        filtered_files = sorted_other_files

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Combine files with improved structure
    with open(output_file, "w", encoding="utf-8") as outfile:
        # Add a comprehensive header
        outfile.write("# CAMB Documentation\n\n")
        outfile.write("---\n\n")

        # Add each file's content
        for file_path in filtered_files:
            file_name = os.path.basename(file_path)
            section_name = os.path.splitext(file_name)[0]
            link_name = 'https://camb.readthedocs.io/en/latest/' + section_name + '.html'

            print(f"  Adding {section_name}...")
            outfile.write(f"## {link_name}\n\n")

            # Add file content
            with open(file_path, "r", encoding="utf-8") as infile:
                content = infile.read()
                outfile.write(content)
                outfile.write("\n\n")
    return True


def convert_CAMBdemo_to_markdown():
    """Convert CAMBdemo.ipynb to markdown using jupytext."""
    print("Converting CAMBdemo.ipynb to markdown...")

    # Path to the notebook
    notebook_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CAMBdemo.ipynb")

    if not os.path.exists(notebook_path):
        print(f"Error: CAMBdemo.ipynb not found at {notebook_path}")
        return None

    # Run jupytext to convert the notebook to markdown
    result = subprocess.run(
        [
            "jupytext",
            "--to", "md",
            "--opt", "notebook_metadata_filter=-all",
            notebook_path
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"Error converting CAMBdemo.ipynb to markdown: {result.stderr}")
        return None

    # Path to the generated markdown file
    md_path = os.path.splitext(notebook_path)[0] + ".md"

    if not os.path.exists(md_path):
        print(f"Error: Generated markdown file not found at {md_path}")
        return None

    print(f"Successfully converted CAMBdemo.ipynb to markdown: {md_path}")
    return md_path


def main():
    args = parse_args()
    try:
        # Get the list of files to exclude
        exclude_files = args.exclude.split(",") if args.exclude else []

        # Build the documentation
        build_dir = build_markdown_docs()

        # Combine the files
        if not combine_markdown_files(build_dir, exclude_files, args.output):
            print("Failed to combine markdown files. Exiting.")
            return 1

        # Verify the file exists and has content
        if not os.path.exists(args.output):
            print(f"ERROR: Failed to generate markdown file at {args.output}")
            return 1

        file_size = os.path.getsize(args.output)
        print(f"Initial markdown file size: {file_size / 1024:.2f} KB")

        if file_size == 0:
            print("ERROR: Generated markdown file is empty")
            return 1

        # Convert CAMBdemo.ipynb to markdown and append it to the combined file
        CAMBdemo_md = convert_CAMBdemo_to_markdown()
        if CAMBdemo_md:
            print(f"Appending CAMBdemo.md to {args.output}...")
            with open(args.output, "a", encoding="utf-8") as outfile:
                outfile.write("## Usage examples from CAMBdemo jupyter notebook\n\n")
                with open(CAMBdemo_md, "r", encoding="utf-8") as infile:
                    content = infile.read()
                    outfile.write(content)

            # Update file size
            file_size = os.path.getsize(args.output)
            print(f"Final markdown file size: {file_size / 1024:.2f} KB")

        print(f"\nSuccess! Documentation has been built and combined into: {args.output}")
        return 0

    except Exception as e:
        print(f"ERROR: An exception occurred during the build process: {e}")
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
