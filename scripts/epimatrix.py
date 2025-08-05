# Import libraries
import argparse
import sys
import os
import glob
from tqdm import tqdm
import dask.dataframe as dd
from dask import delayed, compute
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

class ChrState:
    def __init__(self, chr, firstbase, lastbase, mnemonic, desc, RGB):
        self.chr = chr
        self.firstbase = firstbase
        self.lastbase = lastbase
        self.mnemonic = mnemonic
        self.desc = desc
        self.RGB = RGB

    def __str__(self):
        return f"[{self.chr}:{self.firstbase}-{self.lastbase}] {self.mnemonic} ({self.desc})"

def get_state_color_map():
    # A consistent color palette for each mnemonic
    state_colors = {
        "TssA": "#FF0000", "TssFlnk": "#FF4500", "TssFlnkU": "#FF4500", "TssFlnkD": "#FF4500",
        "Tx": "#008000", "TxWk": "#006400", "EnhG1": "#C2E105", "EnhG2": "#C2E105", "EnhA1": "#FFC34D",
        "EnhA2": "#FFC34D", "EnhWk": "#FFFF00", "ZNF/Rpts": "#66CDAA", "Het": "#8A91D0", "TssBiv": "#CD5C5C",
        "EnhBiv": "#BDB76B", "ReprPC": "#808080", "ReprPCWk": "#C0C0C0", "Quies": "#FFFFFF"
    }
    return state_colors


def MnemonicToDesc(mnemonic):
    # Short description for each mnemonic
    descs = {
        "TssA": "Active TSS", "TssFlnk": "Flanking TSS", "TssFlnkU": "Flanking TSS Upstream",
        "TssFlnkD": "Flanking TSS Downstream", "Tx": "Strong transcription", "TxWk": "Weak transcription",
        "EnhG1": "Genic Enhancer 1", "EnhG2": "Genic Enhancer 2", "EnhA1": "Active Enhancer 1",
        "EnhA2": "Active Enhancer 2", "EnhWk": "Weak Enhancer", "ZNF/Rpts": "ZNF genes & repeats",
        "Het": "Heterochromatin", "TssBiv": "Bivalent/Poised TSS", "EnhBiv": "Bivalent Enhancer",
        "ReprPC": "Repressed PolyComb", "ReprPCWk": "Weak Repressed PolyComb", "Quies": "Quiescent/Low"
    }
    return descs.get(mnemonic, "")


def parse_entries(entry_arg):
    entries = entry_arg.replace(":", "-").split(",")
    return [(e.split("-")[0], int(e.split("-")[1]), int(e.split("-")[2])) for e in entries]


def get_segments_for_window(chrom, start, end, window):
    # Divide the genome into segments of equal size, specified as 'window'
    segments = []
    while start < end:
        segment_end = min(start + window, end)
        segments.append((chrom, start, segment_end))
        start = segment_end
    return segments


@delayed
def process_segment_from_ddf(df, chrom, start, end):
    result = df[(df["chromosome"] == chrom) &
                (df["start"] <= end) &
                (df["end"] >= start)]

    if result.empty:
        return []

    # Determine the majority state for each segment
    chrom_state_counter = result["mnemonic"].value_counts()
    majority_state = chrom_state_counter.idxmax() if not chrom_state_counter.empty else None

    if majority_state:
        desc = MnemonicToDesc(majority_state)
        RGB = result["color"].iloc[0]
        return [ChrState(chrom, start, end, majority_state, desc, RGB)]

    return []


def process_file_with_dask(file_path, segment_list):
    # File processing with dask dataframe
    try:
        ddf = dd.read_csv(
            file_path,
            sep="\t",
            header=None,
            compression='gzip',
            blocksize=None,
            names=["chromosome", "start", "end", "mnemonic", "zero", "dot", "n", "a", "color"],
            usecols=["chromosome", "start", "end", "mnemonic", "color"]
        )
    except Exception as e:
        print(f"Failed to read {file_path}: {e}")
        return []

    try:
        df = ddf.compute()
    except Exception as e:
        print(f"Error computing Dask DataFrame for {file_path}: {e}")
        return []

    state_list = []
    results = []

    for chrom, start, end in tqdm(segment_list, desc=f"Preparing segments from {os.path.basename(file_path)}"):
        results.append(process_segment_from_ddf(df, chrom, start, end))

    # Compute results
    computed_results = list(tqdm(compute(*results), total=len(results), desc="Computing"))

    for result in computed_results:
        if result:
            state_list.extend(result)

    return state_list

# Define sample and state class load functions
def load_sample_classes(classmap_path):
    df = pd.read_csv(classmap_path, sep="\t", header=None, names=["SampleID", "Class"])
    return dict(zip(df["SampleID"], df["Class"]))

def load_state_map(statemap_path):
    df = pd.read_csv(statemap_path, sep="\t", header=None, names=["Mnemonic", "StateClass"])
    return dict(zip(df["Mnemonic"], df["StateClass"]))

# Parse arguments
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", help="Path to the directory containing the raw data", type=str)
    parser.add_argument("entry", help="Chromosome region(s), e.g., 'chr7:140000-150000, chr10:100000-150000'", type=str)
    parser.add_argument("window", help="Window size", type=int)
    parser.add_argument("--classmap", help="Path to sample classification .tsv file", required=True, type=str)
    parser.add_argument("--statemap", help="Path to state classification .tsv file", required=False, type=str)
    parser.add_argument("--output", help="Output file for matrix (.tsv)", default=None)
    parser.add_argument("--plot", help="Output file for state visualization plot (.png)", default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Parse input regions
    parsed_entries = parse_entries(args.entry)
    segment_list = []
    for chrom, start, end in parsed_entries:
        segment_list.extend(get_segments_for_window(chrom, start, end, args.window))

    # Load sample and state classes
    sample_class_map = load_sample_classes(args.classmap)
    state_map = load_state_map(args.statemap) if args.statemap else None

    # Load original color map
    original_state_color_map = get_state_color_map()

    # Map reduced state colors if statemap is provided
    if state_map:
        def get_reduced_state_color_map(state_map, original_state_colors):
            category_colors = {}
            category_to_mnemonics = {}

            # Agrupar Mnemonic por categoría reducida
            for mnemonic, category in state_map.items():
                category_to_mnemonics.setdefault(category, []).append(mnemonic)
            
            for category, mnemonics in category_to_mnemonics.items():
                colors = [original_state_colors.get(m, "#D3D3D3") for m in mnemonics if m in original_state_colors]
                if colors:
                    # Usamos el color del primer Mnemonic que aparece en la categoría
                    category_colors[category] = colors[0]
                else:
                    category_colors[category] = "#D3D3D3"  # Color por defecto

            return category_colors

        reduced_state_color_map = get_reduced_state_color_map(state_map, original_state_color_map)
    else:
        reduced_state_color_map = original_state_color_map

    # Scanning sample files
    gz_files = glob.glob(os.path.join(args.dir, "BSS*_18_CALLS_segments.bed.gz"))

    # Group samples by class
    class_grouped_results = defaultdict(dict)

    for file_path in gz_files:
        sample_id = os.path.basename(file_path).split("_")[0]
        sample_class = sample_class_map.get(sample_id, "Unknown")

        try:
            state_list = process_file_with_dask(file_path, segment_list)
            if state_list:
                class_grouped_results[sample_class][sample_id] = state_list
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    # Write results
    print("Outputting results...\n")
    for class_label, samples in class_grouped_results.items():
        print(f"\n### Class: {class_label} ###\n")
        for sample_id, states in samples.items():
            print(f"--- Sample: {sample_id} ---")
            for state in states:
                print(state)

    # Aggregate state counts per class
    print("\nGenerating state count matrix...\n")
    state_count_matrix = defaultdict(lambda: defaultdict(int))

    for class_label, samples in class_grouped_results.items():
        for sample_id, states in samples.items():
            for state in states:
                # Usa la categoría reducida si está definido
                state_key = state_map.get(state.mnemonic, state.mnemonic) if state_map else state.mnemonic
                state_count_matrix[state_key][class_label] += 1

    # Convert to pandas DataFrame
    matrix_df = pd.DataFrame(state_count_matrix).T.fillna(0).astype(int)
    matrix_df = matrix_df.sort_index(axis=0).sort_index(axis=1)

    print("\n### Chromatin State Count Matrix ###\n")
    print(matrix_df)

    # Optionally write to file
    if args.output:
        matrix_df.to_csv(args.output, sep="\t")
        print(f"\nMatrix saved to: {args.output}")

    # State visualization plot if requested
    if args.plot:
        print(f"\nGenerating state visualization plot at: {args.plot}")

        # Usa el mapa de color reducido si se pasó statemap
        state_color_map = reduced_state_color_map
        default_color = "#D3D3D3"

        segment_label_list = [f"{chrom}:{start}-{end}" for chrom, start, end in segment_list]
        color_matrix = []
        sample_labels = []

        for class_label, samples in class_grouped_results.items():
            for sample_id, states in samples.items():
                row_colors = []
                for chrom, start, end in segment_list:
                    matched = next(
                        (s for s in states if s.chr == chrom and s.firstbase == start and s.lastbase == end),
                        None
                    )
                    # Si statemap está, mapea el mnemonic a categoría
                    if state_map:
                        state = state_map.get(matched.mnemonic, "NA") if matched else "NA"
                    else:
                        state = matched.mnemonic if matched else "NA"
                    
                    color = state_color_map.get(state, default_color)
                    row_colors.append(color)
                color_matrix.append(row_colors)
                sample_labels.append(f"{class_label}:{sample_id}")

        # Safety check: is there data?
        if not color_matrix or not sample_labels:
            print("No samples or color matrix data to plot. Skipping plot.")
            return

        print(f"Total samples to plot: {len(sample_labels)}")
        print(f"Total segments per sample: {len(segment_list)}")

        # Plot
        max_pixels = 65500  # Under the 2^16 limit to be safe
        dpi = 300
        height_per_sample = 0.15
        base_height = len(sample_labels) * height_per_sample + 1.5
        max_height_inches = max_pixels / dpi
        fig_height = min(base_height, max_height_inches)
        fig_width = max(12, len(segment_list) * 0.3)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

        for y, row in enumerate(color_matrix):
            for x, color in enumerate(row):
                ax.add_patch(plt.Rectangle((x, y), 1, 1, color=color, linewidth=0))

        ax.set_xlim(0, len(segment_list))
        ax.set_ylim(0, len(sample_labels))
        ax.set_xticks([i + 0.5 for i in range(len(segment_list))])
        ax.set_yticks([i + 0.5 for i in range(len(sample_labels))])
        ax.set_xticklabels(segment_label_list, rotation=90, fontsize=8)
        ax.set_yticklabels(sample_labels, fontsize=8)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Legend
        used_states = set()
        for row in color_matrix:
            for color in row:
                for state, color_code in state_color_map.items():
                    if color == color_code:
                        used_states.add(state)
        legend_patches = [
            mpatches.Patch(color=color, label=state)
            for state, color in state_color_map.items() if state in used_states
        ]
        ax.legend(handles=legend_patches, bbox_to_anchor=(1.01, 1), loc='upper left', title="States", fontsize=8)

        plt.tight_layout()
        plt.savefig(args.plot, dpi=100, facecolor='white')
        plt.close()
        print(f"Plot saved as: {args.plot}")


if __name__ == "__main__":
    main()