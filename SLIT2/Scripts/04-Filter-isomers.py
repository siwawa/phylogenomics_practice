import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time

metadata_file = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.txt"
fasta_file = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.fasta"
output_fasta = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/SLIT-homologs.fasta"

print("0. 메타데이터 파일에서 Accession과 Organism 정보 매핑 중...")
# 공백/탭 상관없이 테이블을 읽어옵니다.
try:
    meta_df = pd.read_csv(metadata_file, sep='\t')
    # Accession(saccver)을 키(Key)로, organism(Clade_Species)을 값(Value)으로 하는 딕셔너리 생성
    acc_to_organism = dict(zip(meta_df['saccver'], meta_df['organism']))
    print(f"-> 총 {len(acc_to_organism)}개의 고유 매핑 정보를 불러왔습니다.\n")
except Exception as e:
    print(f"메타데이터 파일을 읽는 데 실패했습니다: {e}")
    exit(1)


print("1. FASTA 파일에서 서열, 길이, Accession 추출 중...")
sequences = {}
current_id = ""

with open(fasta_file, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # '>XP_019635574.1 ...' 에서 오직 Accession 번호만 추출
            current_id = line.split()[0][1:]
            sequences[current_id] = ""
        else:
            sequences[current_id] += line

print(acc_to_organism)

# 데이터프레임 생성 (메타데이터 딕셔너리에서 organism 매핑)
data = []
for acc, seq in sequences.items():
    # 딕셔너리에 매핑 정보가 있으면 가져오고, 없으면 예외 처리
    organism_name = acc_to_organism.get(acc, "Unknown_Clade_Species")
    data.append({
        "saccver": acc, 
        "length": len(seq), 
        "sequence": seq,
        "organism": organism_name
    })

df = pd.DataFrame(data)
print(f"-> 총 {len(df)}개의 단백질 서열을 확인했습니다.\n")


print("2. NCBI 서버에서 각 단백질의 Gene ID 추적 중 (약 1~2분 소요)...")
gene_ids = []
for index, acc in enumerate(df['saccver']):
    base_acc = acc.split('.')[0] 
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id={base_acc}"
    
    try:
        response = requests.get(url)
        root = ET.fromstring(response.text)
        
        link = root.find(".//LinkSetDb/Link/Id")
        if link is not None:
            gene_id = link.text
        else:
            print(f"   -> Info: {acc}는 등록된 Gene ID가 없습니다 (독립 서열로 취급).")
            gene_id = f"NoGeneID_{acc}" 
            
    except Exception as e:
        gene_id = f"Error_{acc}"
        
    gene_ids.append(gene_id)
    
    if (index + 1) % 10 == 0:
        print(f"   [{index + 1}/{len(df)}] 매핑 완료...")
    time.sleep(0.4) 

df['GeneID'] = gene_ids


print("\n3. Isoform 제거 및 Paralog 넘버링 (H1, H2...) 적용 중...")
# 1) GeneID로 묶고 가장 긴 서열(Isoform 대표) 1개만 남기기
filtered_df = df.sort_values('length', ascending=False).groupby('GeneID').head(1)

# 2) organism(분류군_학명) 기준으로 그룹화하여 넘버링(H1, H2...) 부여
filtered_df = filtered_df.sort_values(by=['organism', 'length'], ascending=[True, False])
filtered_df['H_num'] = filtered_df.groupby('organism').cumcount() + 1
filtered_df['H_label'] = 'H' + filtered_df['H_num'].astype(str)

print(f"-> 필터링 결과: {len(df)}개 -> {len(filtered_df)}개로 정리되었습니다.\n")


print("4. 최종 FASTA 파일 저장 중...")
with open(output_fasta, "w") as f_out:
    for _, row in filtered_df.iterrows():
        # 새 FASTA 헤더 포맷: >Target_Clade_Binomial_name_H1_geneID original_acc:XP_...
        header = f">{row['organism']}_{row['H_label']}_{row['GeneID']} original_acc:{row['saccver']}"
        f_out.write(header + "\n")
        
        seq = row['sequence']
        for i in range(0, len(seq), 80):
            f_out.write(seq[i:i+80] + "\n")

print(f"🎉 모든 작업이 완료되었습니다! 파일 위치: {output_fasta}")