import os
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
from pathlib import Path


base_dir = Path("/home/lihengtao/data1/prokaryotes/orphan_SelD/final/up_down_15_4/gene_feature_for_paper_20260115")
tsv_dir = base_dir / "tsv"
output_dir = base_dir / "image_pdf"
summary_file = base_dir / "protein_info_summary.xlsx"
colour_file = base_dir / "colour.xlsx"


output_dir.mkdir(parents=True, exist_ok=True)


summary_df = pd.read_excel(summary_file, dtype=str)
colour_df = pd.read_excel(colour_file, dtype=str).set_index("protein")["colour"].to_dict()


for tsv_file in tsv_dir.glob("*.tsv"):
    df = pd.read_csv(tsv_file, sep='\t', dtype=str)
    features = []
    start_list = []
    end_list = []

    for _, row in df.iterrows():
        protein_info = row['protein_info']
        parts = protein_info.split('#')
        if len(parts) < 4:
            continue  

        start = int(parts[1])
        end = int(parts[2])
        frame = int(parts[3])
        strand = 1 if frame > 0 else -1
        label = None
        color = "#ffffff"  

        
        match_found = False
        for col in summary_df.columns:
            if protein_info in summary_df[col].values:
                label = col
                color = colour_df.get(label, "#ffffff")
                match_found = True
                break

        feature_args = {
            "start": start,
            "end": end,
            "strand": strand,
            "color": color
        }
        if match_found:
            feature_args["label"] = label

        features.append(GraphicFeature(**feature_args))
        start_list.append(start)
        end_list.append(end)

    if not features:
        print(f"文件 {tsv_file.name} 无有效蛋白信息，跳过")
        continue

    min_start = min(start_list)
    max_end = max(end_list)
    sequence_length = max_end - min_start + 1800
    first_index = min_start - 900

    
    record = GraphicRecord(sequence_length=sequence_length, first_index=first_index, features=features, ticks_resolution=2000)
    ax, _ = record.plot(figure_width=20, annotate_inline=False)

    
    ax.set_xticks([])  
    ax.set_yticks([])  

    
    ax.figure.canvas.draw()

    # output_path = output_dir / f"{tsv_file.stem}.png"
    # ax.figure.savefig(output_path)
    # print(f"已保存图像：{output_path}")

    
    output_path = output_dir / f"{tsv_file.stem}.pdf"
    ax.figure.savefig(output_path, format='pdf')
    print(f"已保存图像：{output_path}")


