import os
import re
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO

from Bio.Seq import Seq
import warnings
from multiprocessing import Pool
from pathlib import Path
import openpyxl   


def read_genome(genome, genome_fna_dir):
    genome_file = os.path.join(genome_fna_dir, genome + '.fna')

    if not os.path.exists(genome_file):
        raise FileNotFoundError(f"File {genome_file} not found")

    with open(genome_file, 'r', encoding='UTF-8') as fna:
        all_seq = ''
        for line in fna:
            if line.startswith('>'):
                line = line.split(' ')[0] + '\n'
                all_seq += '\n' + line
            else:
                line = line.strip()
                all_seq += line
        all_seq = all_seq.lstrip() + '\n'
        all_list = all_seq.split('\n')
        dict_dbname = {item: all_list[all_list.index(item) + 1]
                       for item in all_list if item.startswith('>')}
        return dict_dbname

def translate(seq):
    seq = seq.strip()
    seq = Seq(seq)
    warnings.filterwarnings("ignore")
    seq = seq.translate()
    return seq

def get_ORF(org_sseq, start, end):
    target_seq = org_sseq[int(start) - 1:int(end)]
    part = str(org_sseq).partition(str(target_seq))
    p1 = part[0][::-1]
    p3 = part[2]

    # 1）寻找上游终止密码子
    for i in range(0, len(p1), 3):
        if i + 2 < len(p1):
            codon = p1[i:i + 3]
            if codon in ['GAT', 'AAT', 'AGT']:
                left_str = p1[:i]
                break

    # 2）如果未找到上游终止密码子
    if 'left_str' in locals():
        pass
    else:
        left_str = p1

    # 3）在该区域内找到最远的起始密码子
    initial_codon_tank = []
    for i in range(0, len(left_str), 3):
        if i + 2 < len(left_str):
            codon = left_str[i:i + 3]
            if codon in ['GTA', 'GTG', 'GTT']:
                initial_codon_tank.append(i + 3)

    if len(initial_codon_tank) == 0:
        target_seq_first3 = target_seq[:3]
        if target_seq_first3 in ['ATG', 'GTG', 'TTG']:
            i_seq = ''
        else:
            if left_str != p1:
                return 'initial_codon not found'
            else:
                left_extra_nucl = len(left_str) % 3
                if left_extra_nucl == 0:
                    i_seq = left_str[::-1]
                elif left_extra_nucl == 1:
                    i_seq = left_str[:-1][::-1]
                elif left_extra_nucl == 2:
                    i_seq = left_str[:-2][::-1]
    else:
        initial_codon = initial_codon_tank[-1]
        i_seq = left_str[:initial_codon][::-1]

    # 4）寻找下游终止密码子
    for i in range(0, len(p3), 3):
        if i + 2 < len(p3):
            codon = p3[i:i + 3]
            if codon in ['TAG', 'TAA', 'TGA']:
                stop_codon = i + 3
                break

    if 'stop_codon' in locals():
        s_seq = p3[:stop_codon]
    else:
        s_seq = p3

    try:
        ORF = i_seq + target_seq + s_seq
        return ORF, int(start) - len(i_seq), int(end) + len(s_seq)
    except:
        return 'ORF not found'

def get_ORF_SEC(org_sseq, all_start):
    target_seq = org_sseq[int(all_start[0]):int(all_start[-1]) + 3]
    p1 = org_sseq[:int(all_start[0])][::-1]
    p3 = org_sseq[int(all_start[-1]) + 3:]

    for i in range(0, len(p1), 3):
        if i + 2 < len(p1):
            codon = p1[i:i + 3]
            if codon in ['GAT', 'AAT', 'AGT']:
                left_str = p1[:i]
                break

    if 'left_str' in locals():
        pass
    else:
        left_str = p1

    initial_codon_tank = []
    for i in range(0, len(left_str), 3):
        if i + 2 < len(left_str):
            codon = p1[i:i + 3]
            if codon in ['GTA', 'GTG', 'GTT']:
                initial_codon_tank.append(i + 3)

    if len(initial_codon_tank) == 0:
        target_seq_first3 = target_seq[:3]
        if target_seq_first3 in ['ATG', 'GTG', 'TTG']:
            i_seq = ''
        else:
            if left_str != p1:
                return 'initial_codon not found'
            else:
                left_extra_nucl = len(left_str) % 3
                if left_extra_nucl == 0:
                    i_seq = left_str[::-1]
                elif left_extra_nucl == 1:
                    i_seq = left_str[:-1][::-1]
                elif left_extra_nucl == 2:
                    i_seq = left_str[:-2][::-1]
    else:
        initial_codon = initial_codon_tank[-1]
        i_seq = left_str[:initial_codon][::-1]

    for i in range(0, len(p3), 3):
        if i + 2 < len(p3):
            codon = p3[i:i + 3]
            if codon in ['TAG', 'TAA', 'TGA']:
                stop_codon = i + 3
                break

    if 'stop_codon' in locals():
        s_seq = p3[:stop_codon]
    else:
        s_seq = p3

    try:
        ORF = i_seq + target_seq + s_seq
        return ORF, int(all_start[0]) - len(i_seq), int(all_start[-1]) + 3 + len(s_seq)
    except:
        return 'ORF not found'

def process_subdir(args):
    subdir, tblastn_dir, genome_fna_dir = args
    result_filedir = os.path.join(tblastn_dir, subdir)

   
    if not os.path.isdir(result_filedir):
        # print(f"Directory {result_filedir} does not exist. Skipping.")
        return

    genomes = os.listdir(result_filedir)

    
    output_subdir = f'{subdir}_ORF'
    output_dir = os.path.join(tblastn_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)

    for geno in genomes:
        if not geno.endswith('_tblastn.tsv'):
            continue
        geno_id = geno.replace('_tblastn.tsv', '')
        # print(f"Processing {geno_id} in {subdir}")
        tblastn_output_path = os.path.join(result_filedir, geno)
        tblastn_output = pd.read_csv(
            tblastn_output_path,
            sep='\t',
            names=[
                'stitle', 'qseqid', 'sseqid', 'qstart', 'qend', 'sstart',
                'send', 'evalue', 'qcovhsp', 'qcovs', 'sframe', 'sseq',
                'qseq', 'pident'
            ]
        )

        ORF_found = False
        result_content = []

        for index, row in tblastn_output.iterrows():
            qseqid = row.iloc[1]
            sstart = row.iloc[5]
            send = row.iloc[6]
            sseq = row.iloc[11]
            qseq = row.iloc[12]

            ori_sseqid = re.compile('[A-Z]*[0-9]*\.[0-9]').findall(row['sseqid'])
            sseqid = '>' + ''.join(ori_sseqid)
            genome_dict = read_genome(geno_id, genome_fna_dir)

            if sseqid not in genome_dict:
                # print(f"{sseqid} not found in genome {geno_id}")
                continue

            if subdir == "SelD_filterbyevalue":
                if int(sstart) < int(send):
                    org_sseq = genome_dict[sseqid]
                    start = int(sstart)
                    end = int(send)
                else:
                    i_seq = genome_dict[sseqid][::-1]
                    org_sseq = Seq(i_seq).reverse_complement()[::-1]
                    length = int(sstart) - int(send)
                    start = len(org_sseq) - int(sstart) + 1
                    end = length + start

                qseq_loca_U = [idx for (idx, value) in enumerate(qseq) if value == 'U']
                sseq_qseq_loca = [sseq[i] for i in qseq_loca_U]

                if '*' in sseq_qseq_loca or 'C' in sseq_qseq_loca:
                    dashcount = [sseq[:i].count('-') for i in qseq_loca_U]
                    real_CU_loca = [qseq_loca_U[i] - dashcount[i] for i in range(len(qseq_loca_U))]
                    ref_U = []
                    for i in range(len(sseq_qseq_loca)):
                        if sseq_qseq_loca[i] == '*' or sseq_qseq_loca[i] == 'C':
                            ref_U.append(real_CU_loca[i])

                    all_start = [int(start) + 2 + i * 3 - 3 for i in ref_U]
                    # print(all_start)
                    if all_start:
                        ORF_result = get_ORF_SEC(org_sseq, all_start)
                    else:
                        ORF_result = 'ORF not found'
                else:
                    # print('no all_start')
                    ORF_result = 'ORF not found'

            else:
                if int(sstart) < int(send):
                    org_sseq = genome_dict[sseqid]
                    start = int(sstart)
                    end = int(send)
                else:
                    i_seq = genome_dict[sseqid][::-1]
                    org_sseq = Seq(i_seq).reverse_complement()[::-1]
                    length = int(sstart) - int(send)
                    start = len(org_sseq) - int(sstart) + 1
                    end = length + start

                ORF_result = get_ORF(org_sseq, start, end)

            if ORF_result == 'initial_codon not found' or ORF_result == 'ORF not found':
                # print(f"ORF not found for {geno_id}, {sseqid}")
                continue
            else:
                ORF_seq, sstart_new, send_new = ORF_result
                ORF_nr = str(translate(ORF_seq))

                if len(ORF_nr) == 0:
                    ORF_nr = ""
                elif ORF_nr[-1] == '*':
                    ORF_nr = ORF_nr[::-1].replace('*', '', 1)[::-1]
                else:
                    ORF_nr = ORF_nr

                ORF_nr = ORF_nr.replace("*", "U")

                # 使用正则表达式替换 CxxK 为 UxxK
                ORF_nr = re.sub(r'C([A-Z]{2})K', r'U\1K', ORF_nr)

                ORF_found = True
                result_content.append(
                    '\t'.join(map(str, row.tolist())) + f'\t{sseqid}\t{ORF_nr}\t{sstart_new}\t{send_new}\n')

        if ORF_found:
            result_path = os.path.join(output_dir, f'{geno_id}.tsv')
            with open(result_path, 'w') as f:
                f.writelines(result_content)

        # print(f"Completed processing {geno_id} in {subdir}.")

    # print(f"Completed processing subdir {subdir}.")

def calculate_overlap(row1, row2):
    """计算两个片段在 query 上的重叠率并返回 (rate1, rate2)。"""
    start1, end1 = row1['qstart'], row1['qend']
    start2, end2 = row2['qstart'], row2['qend']
    
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start > overlap_end:
        return 0, 0  # 无重叠
    
    overlap_length = overlap_end - overlap_start + 1
    len1 = end1 - start1 + 1
    len2 = end2 - start2 + 1
    
    rate1 = overlap_length / len1 if len1 != 0 else 0
    rate2 = overlap_length / len2 if len2 != 0 else 0
    
    return rate1, rate2


def parse_and_group(lines):
    """将rpsblast结果按qaccver和query位点重叠进行复杂分组"""
    
    columns = ['stitle', 'qaccver', 'qstart', 'qend', 'pident', 'qlen', 'length', 'evalue', 'bitscore']
    data = []
    

    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) != 9:
            continue
        evalue = float(parts[7])
        if evalue <= 1e-10:
            data.append({
                'stitle': parts[0],
                'qaccver': parts[1],
                'qstart': int(parts[2]),
                'qend': int(parts[3]),
                'pident': float(parts[4]),
                'qlen': int(parts[5]),
                'length': int(parts[6]),
                'evalue': evalue,
                'bitscore': float(parts[8]),
                'raw': parts  # 保留原始列表内容
            })

    df = pd.DataFrame(data)
    if df.empty:
        return {}

    grouped_lines = {}

    for qaccver, group_df in df.groupby('qaccver'):
        group_df = group_df.sort_values(by='evalue')  # 从最小evalue开始
        remaining = group_df.copy()

        while not remaining.empty:
            # 选取evalue最小的为种子
            seed_row = remaining.iloc[0]
            seed_start, seed_end = seed_row['qstart'], seed_row['qend']
            seed_id = f"{qaccver}_{seed_start}_{seed_end}"

            # 与种子比对，找出满足重叠条件的行
            selected_indices = []
            for idx, row in remaining.iterrows():
                rate1, rate2 = calculate_overlap(seed_row, row)
                if (rate1 > 0.7 or rate2 > 0.7) or (rate1 > 0.5 and rate2 > 0.5):
                    selected_indices.append(idx)

            
            selected_rows = remaining.loc[selected_indices]
            grouped_lines[seed_id] = [row['raw'] for _, row in selected_rows.iterrows()]

            
            remaining = remaining.drop(index=selected_indices)

    return grouped_lines


def run_rpsblast(args):
    fasta_file, result_subdir = args
    accession = fasta_file.stem  
    output_file = result_subdir / f"{accession}.out"


    subprocess.run([
        "rpsblast",
        "-query", str(fasta_file),
        "-db", "/home/lihengtao/data/prokaryotes/db/Cdd",
        "-out", str(output_file),
        "-outfmt", "6 stitle qaccver qstart qend pident qlen length evalue bitscore"
    ])

    # print(f"Processed {fasta_file} -> {output_file}")


def step1_tblastn(query_dir, genome_dir, output_dir, base_dir, classification):

    task_file = os.path.join(base_dir, f"tasks_{classification}.txt")
    with open(task_file, 'w') as tf:

        for query_file in Path(query_dir).glob("*.fasta"):
            query_name = query_file.stem  # 获取蛋白质名称（去掉后缀 .fasta）
            query_output_dir = output_dir / query_name


            query_output_dir.mkdir(parents=True, exist_ok=True)


            for genome_subdir in Path(genome_dir).iterdir():
                if genome_subdir.is_dir():
                    genome_dir_name = genome_subdir.name  
                    genome_db = genome_subdir / f"{genome_dir_name}_db"


                    tblastn_command = (
                        f"tblastn -query {query_file} -db {genome_db} "
                        f"-out {query_output_dir}/{genome_dir_name}_tblastn.tsv "
                        f"-outfmt '6 stitle qseqid sseqid qstart qend sstart send evalue qcovhsp qcovs sframe sseq qseq pident'\n"
                    )
                    tf.write(tblastn_command)


    with open(task_file, 'r') as tf:
        subprocess.run(['parallel', '-j', '80'], stdin=tf)

    print("步骤一：tblastn分析完成。")

def step2_filter_evalue(tblastn_dir, round_number):

    if round_number == 1:
        evalue_threshold = 5
    else:
        evalue_threshold = 5


    sub_dirs = [d for d in os.listdir(tblastn_dir) if os.path.isdir(os.path.join(tblastn_dir, d))]


    def process_tsv(tsv_file, output_file):
        try:

            if os.stat(tsv_file).st_size == 0:
                # print(f"{tsv_file} 文件为空，跳过处理。")
                return


            df = pd.read_csv(tsv_file, sep='\t', header=None)

            df.columns = ['stitle', 'qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send',
                          'evalue', 'qcovhsp', 'qcovs', 'sframe', 'sseq', 'qseq', 'pident']


            df_filtered = df[df['evalue'] <= evalue_threshold]


            df_final = df_filtered.loc[df_filtered.groupby('sseqid')['evalue'].idxmin()]


            df_final.to_csv(output_file, sep='\t', header=False, index=False)

        except pd.errors.EmptyDataError:
            print(f"警告: {tsv_file} 为空或格式不正确，跳过处理。")
        except Exception as e:
            print(f"处理文件 {tsv_file} 时发生错误: {e}")


    for sub_dir in sub_dirs:
        sub_dir_path = os.path.join(tblastn_dir, sub_dir)
        if not os.path.isdir(sub_dir_path):
            continue


        output_dir = os.path.join(tblastn_dir, f"{sub_dir}_filterbyevalue")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)


        for file in os.listdir(sub_dir_path):
            if file.endswith('_tblastn.tsv'):
                file_path = os.path.join(sub_dir_path, file)
                output_file_path = os.path.join(output_dir, file)


                process_tsv(file_path, output_file_path)

    print("步骤二：根据E值过滤完成。")


def step3_find_ORF(tblastn_dir, genome_fna_dir):
    subdirs = [
        'OvsA_strict_filterbyevalue', 'OvsM_filterbyevalue', 'SelA_filterbyevalue',
        'SelB_filterbyevalue', 'SelD_filterbyevalue', 'SenA_strict_filterbyevalue',
        'SenB_filterbyevalue', 'YbbB_filterbyevalue', 'YqeB_filterbyevalue',
        'YqeC_filterbyevalue', 'OvsA_lenient_filterbyevalue', 'SenA_lenient_filterbyevalue'
    ]


    existing_subdirs = [subdir for subdir in subdirs if os.path.isdir(os.path.join(tblastn_dir, subdir))]


    with Pool(processes=12) as pool:
        args = [(subdir, tblastn_dir, genome_fna_dir) for subdir in existing_subdirs]
        pool.map(process_subdir, args)

    print("步骤三：ORF查找完成。")


def step4_write_ORF_to_fasta(tblastn_dir):
    subdirs = [
        "OvsA_strict_filterbyevalue_ORF", "OvsM_filterbyevalue_ORF", "SelA_filterbyevalue_ORF",
        "SelB_filterbyevalue_ORF", "SelD_filterbyevalue_ORF", "SenA_strict_filterbyevalue_ORF",
        "SenB_filterbyevalue_ORF", "YbbB_filterbyevalue_ORF", "YqeB_filterbyevalue_ORF",
        "YqeC_filterbyevalue_ORF", "OvsA_lenient_filterbyevalue_ORF", "SenA_lenient_filterbyevalue_ORF"
    ]

    for subdir in subdirs:
        subdir_path = os.path.join(tblastn_dir, subdir)

        if not os.path.isdir(subdir_path):
            continue

        hitted_dir = os.path.join(subdir_path, 'hitted')
        os.makedirs(hitted_dir, exist_ok=True)

        for filename in os.listdir(subdir_path):
            if filename.endswith(".tsv"):
                accession = filename.replace(".tsv", "")
                fasta_file = os.path.join(hitted_dir, f"{accession}.fasta")

                with open(os.path.join(subdir_path, filename), 'r') as tsv_file, open(fasta_file, 'w') as fasta_out:
                    for line in tsv_file:
                        columns = line.strip().split("\t")

                        if len(columns) < 17:
                            continue

                        seq_id = columns[-4].replace(">", "") 
                        orf_seq = columns[-3]  
                        start_pos = columns[-2]  
                        end_pos = columns[-1]  
                        frame_info = columns[10]  

                        header = f">{seq_id}_{start_pos}_{end_pos}_{frame_info}"
                        sequence = orf_seq

                        fasta_out.write(f"{header}\n{sequence}\n")

    print("步骤四：ORF写入fasta文件完成。")


def step5_generate_ORF_excel(tblastn_dir, output_excel_path):
    
    prefixes = ['SelA', 'SelB', 'SelD', 'YbbB', 'YqeB', 'YqeC', 'SenA_strict', 'SenB', 'OvsA_strict', 'OvsM',
                'SenA_lenient', 'OvsA_lenient']

    
    results = []

    
    for prefix in prefixes:
        dir_path = os.path.join(tblastn_dir, f"{prefix}_filterbyevalue_ORF")

        if not os.path.exists(dir_path):
            continue

        
        tsv_files = [f for f in os.listdir(dir_path) if f.endswith('.tsv')]
        num_files = len(tsv_files)

        total_lines = 0
        for tsv_file in tsv_files:
            file_path = os.path.join(dir_path, tsv_file)
            with open(file_path, 'r') as f:
                total_lines += sum(1 for _ in f)

        
        results.append([prefix, num_files, total_lines])

    
    df = pd.DataFrame(results, columns=['Directory', 'TSV File Count', 'Total Lines'])

    
    df.to_excel(output_excel_path, index=False)

    print(f"步骤五：ORF统计结果已写入 {output_excel_path}")


def step5_rpsblast(tblastn_dir, rpsblast_dir):
    
    subdirs = ['SelA', 'SelB', 'SelD', 'YbbB', 'YqeB', 'YqeC', 'SenA_strict', 'SenB', 'OvsA_strict', 'OvsM',
               'SenA_lenient', 'OvsA_lenient']

    
    results_dir = os.path.join(rpsblast_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    for subdir in subdirs:
        os.makedirs(os.path.join(results_dir, subdir), exist_ok=True)

    
    for subdir in subdirs:
        query_dir = os.path.join(tblastn_dir, f"{subdir}_filterbyevalue_ORF/hitted")
        result_subdir = Path(results_dir) / subdir

        if not os.path.isdir(query_dir):
            continue

        
        fasta_files = list(Path(query_dir).glob("*.fasta"))

        with Pool(processes=80) as pool:
            args = [(fasta_file, result_subdir) for fasta_file in fasta_files]
            pool.map(run_rpsblast, args)

    print("步骤六：rpsblast分析完成！")


def step6_parse_rpsblast_results(rpsblast_results_dir, tblastn_dir, output_excel):
    # 每种蛋白质对应的domain号
    domains = {
        "SelD": ["cd02195", "PRK14105", "TIGR00476", "PRK00943", "COG0709"],
        "SelA": ["TIGR00474", "pfam03841", "COG1921"],
        "SelB": ["PRK10512", "TIGR00475", "COG3276", "cd04171"],
        "YbbB": ["PRK11784", "COG2603", "TIGR03167"],
        "YqeB": ["TIGR03309"],
        "YqeC": ["TIGR03172", "pfam19842"],
        "SenB": ["TIGR04348"],
        "SenA_lenient": ["NF041186"],
        "OvsA_lenient": ["NF041186"],
        "SenA_strict": ["NF041186"],
        "OvsA_strict": ["NF041186"],
        "OvsM": ["TIGR04345"]
    }

    
    with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
        
        for protein, valid_domains in domains.items():
            protein_dir = os.path.join(rpsblast_results_dir, protein)
            if not os.path.isdir(protein_dir):
                continue

            result_dict = {}

            
            for file_name in os.listdir(protein_dir):
                if file_name.endswith(".out"):
                    accession = file_name.replace(".out", "")  # 提取accession
                    file_path = os.path.join(protein_dir, file_name)

                    
                    with open(file_path, "r") as file:
                        lines = file.readlines()
                    
                    grouped_lines = parse_and_group(lines)

                    # 检查domain号并记录qaccver
                    qaccver_set = set()
                    if protein == "OvsM":
                        invalid_domain = "TIGR04344"
                        # delete_qaccver = set()
                        for qaccver, entries in grouped_lines.items():
                            # 取最小的三行
                            top_entries = sorted(entries, key=lambda x: float(x[7]))[:3]
                            domain_numbers = {entry[0].split(',')[0].strip() for entry in top_entries}
                            domain_text = top_entries[0][0].split(',')[0].strip()
                            if invalid_domain in domain_numbers:
                                continue
                            elif domain_text in valid_domains:
                                qaccver_set.add(qaccver)

                    else:
                        for qaccver, entries in grouped_lines.items():
                            top_entry = min(entries, key=lambda x: float(x[7]))
                            domain_text = top_entry[0].split(',')[0].strip()
                            if domain_text in valid_domains:
                                qaccver_set.add(qaccver)
                        if protein == "SenA_strict":
                            # print(accession)
                            final_qaccver_set = set()
                            for acc in qaccver_set:
                                acc_ture = acc.rsplit('_', 2)[0]
                                fasta_file_path = os.path.join(tblastn_dir,
                                                               f"{protein}_filterbyevalue_ORF/hitted/{accession}.fasta")
                                # print(fasta_file_path)
                                if os.path.isfile(fasta_file_path):

                                    print(fasta_file_path)
                                    found = False
                                    for record in SeqIO.parse(fasta_file_path, "fasta"):
                                        if record.id == acc_ture:
                                            print(acc)
                                            sequence = str(record.seq)
                                            found = True
                                            print(sequence)
                                            break
                                    if found and "MH" in sequence and any(
                                            motif in sequence for motif in ["NFF", "NFY", "NYY", "NYF"]):
                                        print("!!!!!!")
                                        final_qaccver_set.add(acc)
                            qaccver_set = final_qaccver_set
                        elif protein == "OvsA_strict":
                            final_qaccver_set = set()
                            for acc in qaccver_set:
                                acc_ture = acc.rsplit('_', 2)[0]
                                fasta_file_path = os.path.join(tblastn_dir,
                                                               f"{protein}_filterbyevalue_ORF/hitted/{accession}.fasta")
                                if os.path.isfile(fasta_file_path):
                                    found = False
                                    for record in SeqIO.parse(fasta_file_path, "fasta"):
                                        if record.id == acc_ture:
                                            sequence = str(record.seq)
                                            found = True
                                            break
                                    if found and "MH" in sequence and "QAY" in sequence:
                                        final_qaccver_set.add(acc)
                            qaccver_set = final_qaccver_set

                    
                    if accession in result_dict:
                        result_dict[accession].update(qaccver_set)
                    else:
                        result_dict[accession] = qaccver_set

            
            data = []
            for accession, qaccver_set in result_dict.items():
                
                row = [accession] + list(qaccver_set)
                max_columns = max([len(x) for x in result_dict.values()]) + 2
                row += [''] * (max_columns - len(row))  
                row.append(len(qaccver_set))  
                data.append(row)

            
            if data:
                df = pd.DataFrame(data)
                df.to_excel(writer, sheet_name=protein, header=False, index=False)

    print(f"步骤七：rpsblast结果已写入 {output_excel}")


def step7_extract_sequences(tblastn_dir, rpsblast_results_excel, output_dir):
    
    xls = pd.ExcelFile(rpsblast_results_excel)

    
    for sheet_name in xls.sheet_names:
        
        protein_dir = os.path.join(output_dir, sheet_name)
        os.makedirs(protein_dir, exist_ok=True)

        
        df = pd.read_excel(xls, sheet_name=sheet_name, header=None)

        
        for index, row in df.iterrows():
            accession = row[0]  
            seq_ids = [s.rsplit('_', 2)[0] for s in row.iloc[1:-1].dropna()]
            seq_count = row.iloc[-1]  

            
            if not seq_ids:
                # print(f"Skipping {accession} because there are no seq_ids.")
                continue

            
            fasta_file_path = os.path.join(tblastn_dir, f'{sheet_name}_filterbyevalue_ORF/hitted/{accession}.fasta')

            
            if os.path.exists(fasta_file_path):
                with open(fasta_file_path, 'r') as fasta_file:
                    fasta_lines = fasta_file.readlines()

                
                output_fasta_path = os.path.join(protein_dir, f'{accession}.fasta')
                with open(output_fasta_path, 'w') as output_file:
                    
                    write_seq = False
                    for line in fasta_lines:
                        if line.startswith('>'):
                            
                            seq_id = line[1:].strip()  
                            if seq_id in seq_ids:
                                output_file.write(line)  
                                write_seq = True  
                            else:
                                write_seq = False  
                        elif write_seq:
                            output_file.write(line)  

    print("步骤八：序列提取完成。")


def step8_combine_fasta_files(output_dir):
    
    for sub_dir in os.listdir(output_dir):
        sub_dir_path = os.path.join(output_dir, sub_dir)

        if os.path.isdir(sub_dir_path):
            combined_file_path = os.path.join(sub_dir_path, 'combined.fasta')

            
            with open(combined_file_path, 'w') as combined_file:
                for fasta_file in os.listdir(sub_dir_path):
                    if fasta_file.endswith('.fasta') and fasta_file != 'combined.fasta':
                        fasta_file_path = os.path.join(sub_dir_path, fasta_file)

                        
                        for record in SeqIO.parse(fasta_file_path, 'fasta'):
                            SeqIO.write(record, combined_file, 'fasta')

    print("步骤九：fasta文件合并完成。")


def step9_run_cdhit(output_dir, next_query_dir):
    
    for protein_name in os.listdir(output_dir):
        subdir_path = os.path.join(output_dir, protein_name)

        
        if os.path.isdir(subdir_path):
            fasta_file = os.path.join(subdir_path, 'combined.fasta')

            if os.path.exists(fasta_file):
                output_file = os.path.join(next_query_dir, f'{protein_name}.fasta')

                
                cmd = [
                    'cd-hit',
                    '-i', fasta_file,
                    '-o', output_file,
                    '-c', '0.4',
                    '-n', '2',
                    '-g', '1'
                ]
                subprocess.run(cmd)

    print("步骤十：cd-hit操作完成。")



def step10_copy_to_final(output_dir, final_dir):
    subdirs = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]

    new_accessions = {}

    for subdir in subdirs:
        source_subdir = os.path.join(output_dir, subdir)
        target_subdir = os.path.join(final_dir, subdir)

        
        os.makedirs(target_subdir, exist_ok=True)

        
        fasta_files = [f for f in os.listdir(source_subdir) if f.endswith('.fasta') and f != 'combined.fasta']

        
        new_accessions[subdir] = []

        
        for file in fasta_files:
            accession = file.replace('.fasta', '')
            source_file = os.path.join(source_subdir, file)
            target_file = os.path.join(target_subdir, file)

            
            if not os.path.exists(target_file):
                new_accessions[subdir].append(accession)
                shutil.copy(source_file, target_subdir)
                # print(f"复制文件 {file} 到 {target_subdir}")

    print("文件复制完成！")
    return new_accessions


def write_accessions_to_excel(accessions_excel, new_accessions, round_number):
    
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in new_accessions.items()]))

    sheet_name = f'Iteration_{round_number}'

    if not os.path.exists(accessions_excel):
        
        with pd.ExcelWriter(accessions_excel, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        
        with pd.ExcelWriter(accessions_excel, engine='openpyxl', mode='a') as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)



def main():
    
    all_bacteria_dir = '/home/lihengtao/data1/prokaryotes/all_bacteria_20250414'
    genome_root_dir = '/home/lihengtao/data/prokaryotes/original/genomes_new/genomes_lht'

    
    bacterial_classifications = [name for name in os.listdir(all_bacteria_dir)
                                 if os.path.isdir(os.path.join(all_bacteria_dir, name))]

    
    for classification in bacterial_classifications:
        print(f"Processing bacterial classification: {classification}")

        base_dir = os.path.join(all_bacteria_dir, classification)
        genome_dir = os.path.join(genome_root_dir, classification)
        genome_fna_dir = os.path.join(genome_root_dir, classification + '_fna')

        final_dir = os.path.join(base_dir, 'final')
        query_base = os.path.join(base_dir, 'query_')

        
        round_number = 1

        
        subdirs = ['SelA', 'SelB', 'SelD', 'YbbB', 'YqeB', 'YqeC', 'SenA_strict', 'SenB',
                   'OvsA_strict', 'OvsM', 'SenA_lenient', 'OvsA_lenient']

        while True:
            print(f"开始第 {round_number} 轮迭代...")

            
            current_query_dir = f"{query_base}{round_number}"
            tblastn_dir = os.path.join(base_dir, f'tblastn_{round_number}')
            output_dir = os.path.join(base_dir, f'seq_{round_number}')
            rpsblast_dir = os.path.join(base_dir, f'rpsblast_{round_number}')
            rpsblast_results_dir = os.path.join(rpsblast_dir, 'results')
            output_excel_path = os.path.join(tblastn_dir, f'ORF_count_{round_number}.xlsx')
            rpsblast_results_excel = os.path.join(rpsblast_results_dir, 'rpsblast_results.xlsx')
            next_query_dir = f"{query_base}{round_number + 1}"

            
            os.makedirs(current_query_dir, exist_ok=True)
            os.makedirs(tblastn_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)
            os.makedirs(rpsblast_results_dir, exist_ok=True)
            os.makedirs(next_query_dir, exist_ok=True)
            os.makedirs(final_dir, exist_ok=True)

            # 步骤一：运行 tblastn
            step1_tblastn(
                query_dir=Path(current_query_dir),
                genome_dir=genome_dir,
                output_dir=Path(tblastn_dir),
                base_dir=base_dir,
                classification=classification
            )

            # 步骤二：根据 E 值过滤 tblastn 结果
            step2_filter_evalue(tblastn_dir, round_number)

            # 步骤三：寻找 ORF
            step3_find_ORF(tblastn_dir, genome_fna_dir)

            # 步骤四：将 ORF 写入 fasta 文件
            step4_write_ORF_to_fasta(tblastn_dir)

            # 步骤五：生成 ORF 统计 Excel 文件
            step5_generate_ORF_excel(tblastn_dir, output_excel_path)

            # 步骤六：运行 rpsblast
            step5_rpsblast(tblastn_dir, rpsblast_dir)

            # 步骤七：解析 rpsblast 结果并写入 Excel
            step6_parse_rpsblast_results(rpsblast_results_dir, tblastn_dir, rpsblast_results_excel)

            # 步骤八：提取具有 domain 的蛋白序列
            step7_extract_sequences(tblastn_dir, rpsblast_results_excel, output_dir)

            # 步骤九：合并 fasta 文件
            step8_combine_fasta_files(output_dir)

            # 步骤十：运行 cd-hit
            step9_run_cdhit(output_dir, next_query_dir)

            # 步骤十一：复制文件到最终目录并检查是否需要继续迭代
            new_accessions = step10_copy_to_final(output_dir, final_dir)

            # 打印每轮新加入的 accessions
            print(f"第 {round_number} 轮迭代，新加入的 accessions:")
            for protein, accs in new_accessions.items():
                print(f"{protein}: {accs}")

            # 将新加入的 accessions 写入 Excel 文件
            accessions_excel = os.path.join(final_dir, 'accessions.xlsx')
            write_accessions_to_excel(accessions_excel, new_accessions, round_number)

            # 判断是否有新的 accession 加入
            if all(len(acc_list) == 0 for acc_list in new_accessions.values()):
                print("所有蛋白均无新加入的 accession，停止迭代。")
                break
            else:
                print("有新的 accession 加入，继续迭代。")
                round_number += 1  # 增加轮次

        print(f"所有操作完成。完成细菌分类：{classification}")

if __name__ == '__main__':
    main()
