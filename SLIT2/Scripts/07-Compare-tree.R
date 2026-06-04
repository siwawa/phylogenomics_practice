
pdf("/rna/liha/phylogenomics_practice/SLIT2/Tree/Raw-alignment/Cladogram.pdf", width = 15, height = 15) 

raw_tree <- read.tree("/rna/liha/phylogenomics_practice/SLIT2/Tree/Raw-alignment/SLIT.treefile") 

# 2. "Branchiostoma"가 포함된 모든 서열 이름(Tip label) 자동 추출
outgroup_tips <- raw_tree$tip.label[grep("Branchiostoma", raw_tree$tip.label, ignore.case = TRUE)]

# 추출된 이름 확인 (제대로 찾아졌는지 콘솔에서 확인하세요!)
print(outgroup_tips) 

# 3. 트리 Rooting (뿌리 내리기)
# ★ 주의: edgelabel = TRUE 를 반드시 넣어야 Bootstrap 값이 엉뚱한 노드로 밀리지 않습니다!
rooted_tree <- root(raw_tree, outgroup = outgroup_tips, resolve.root = TRUE, edgelabel = TRUE)

# 4. Rooting 된 트리 시각화
plot(rooted_tree, use.edge.length = FALSE, main = "SLIT phylogeny", cex = 0.6)

# 5. Bootstrap 값 덧그리기 (이전과 동일하게 90 이상만 빨간색으로 표시)
bs_values <- as.numeric(rooted_tree$node.label)
nodelabels(bs_values, cex = 0.5, frame = "none", col = "red", adj = c(1.2, -0.5))

# add.scale.bar(x = 0.1, y = 5, cex = 0.7, lwd = 2)


dev.off()