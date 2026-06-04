import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time
import re

fasta_file = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.fasta"
output_fasta = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-filtered-longest-isoforms-renamed.fasta"

print("1. FASTA 파일에서 서열, 길이, 종(Species) 학명 추출 중...")
sequences = {}
species_dict = {}
current_id = ""

with open(fasta_file, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # Accession 번호 추출
            current_id = line.split()[0][1:]
            sequences[current_id] = ""
            
            # 정규식을 이용해 대괄호 [ ] 안의 학명(Binomial name) 추출
            match = re.search(r'\[(.*?)\]', line)
            if match:
                # 빈칸을 언더바(_)로 변경 (예: Aplidium turbinatum -> Aplidium_turbinatum)
                species_dict[current_id] = match.group(1).replace(' ', '_')
            else:
                species_dict[current_id] = "Unknown_species"
        else:
            sequences[current_id] += line

# 데이터프레임 생성 (species 컬럼 추가)
data = [{
    "saccver": acc, 
    "length": len(seq), 
    "sequence": seq,
    "species": species_dict[acc]
} for acc, seq in sequences.items()]
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
            # 에러 대신 정보성 메시지로 출력
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

# 2) 종(Species)별로 그룹화하여 넘버링(H1, H2...) 부여
# (길이가 긴 서열이 우선적으로 H1을 받도록 정렬)
filtered_df = filtered_df.sort_values(by=['species', 'length'], ascending=[True, False])
filtered_df['H_num'] = filtered_df.groupby('species').cumcount() + 1
filtered_df['H_label'] = 'H' + filtered_df['H_num'].astype(str)

print(f"-> 필터링 결과: {len(df)}개 -> {len(filtered_df)}개로 정리되었습니다.\n")


print("4. 최종 FASTA 파일 저장 중...")
with open(output_fasta, "w") as f_out:
    for _, row in filtered_df.iterrows():
        # 새 FASTA 헤더 포맷: >Binomial_name_H1_geneID original_acc:XP_...
        # 나중에 문제가 생겼을 때 원본을 추적할 수 있도록 원래 Accession 번호는 띄어쓰고 뒤에 달아둡니다.
        header = f">{row['species']}_{row['H_label']}_{row['GeneID']} original_acc:{row['saccver']}"
        f_out.write(header + "\n")
        
        seq = row['sequence']
        for i in range(0, len(seq), 80):
            f_out.write(seq[i:i+80] + "\n")

print(f"🎉 모든 작업이 완료되었습니다! 파일 위치: {output_fasta}")