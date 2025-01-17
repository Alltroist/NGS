####################################
#
#  �������� ������� 3. RNA-Seq
#
#       ������� ������
#
####################################


library(edgeR)

 
######------------------------------------------
#

#                    ����� 1

# ��������� ������ �� ����������� �� ��� ������ edgeR)
# ����� ���� � ������������ �/��� �����������
# ����������� ��������� ���������������� p-value < 0.05,
# ����������� ������� ������ ���� �� ����� ��� � ���
# ����), ���������� ��������� ����� ������� ��������� ��
# ��������. ����������� ������ ~ tissue + age

#
######------------------------------------------


# � count.out �� ������������ �������� �������������
# ������� ������� �� �������

col_names = c("B14.5", "B15.5", "B17.5", "B20", "B34",
              "C14.5", "C15.5", "C17.5", "C20", "C34")

counts = read.table(file = 'count.out', stringsAsFactors = FALSE,
                    col.names = c("gene",col_names))

# ������ ����, ������� �� ������������� �� ���� ��������
counts = counts[rowSums(counts[, 2:11]) >= 1, ] 

# � ��������� ������������ ���������� � ���� ����� � ��������
tissue = as.factor(substring(col_names, 1, 1))
age = as.numeric(substring(col_names, 2))

# �������� � ��������� � ������� ������ edgeR
# ������� �� �������� � ���������� ��������
edger = DGEList(counts = as.matrix(counts[,2:11]), group = tissue)
edger = calcNormFactors(edger, method = 'RLE') 
norm_counts = cpm(edger) 

# ����� ������������ ������ ~tissue + age
design = model.matrix(~ tissue + age)
edger = estimateGLMCommonDisp(edger,design)
edger = estimateGLMTrendedDisp(edger,design)
edger = estimateGLMTagwiseDisp(edger,design)
strict.disp = pmax(edger$tagwise.dispersion, edger$trended.dispersion, edger$common.dispersion)
plotBCV(edger)

#������ �������� � ���������� � �������� GLM ������ 
glm = glmFit(edger, design, dispersion=strict.disp)
lrt_tissue = glmLRT(glm, 2)
lrt_age = glmLRT(glm, 3)

# ����������� ������ � ������� topTags
filt_age = topTags(lrt_age, n = Inf, adjust.method = 'BH', sort.by = 'none')
filt_tissue = topTags(lrt_tissue, n = Inf, adjust.method = 'BH', sort.by = 'none')

# ������ ������ �������� �� ���� ����� ���� 
# p-val ������ ���� < 0.05 � ����������� ������� 
# ������ ���� �� ����� ��� � ��� ����
sign_gene_tissue = (filt_tissue$table$logFC >= 1 | filt_tissue$table$logFC <= -1) & (filt_tissue$table$PValue < 0.05)

#�� ����������: 73
sum(sign_gene_tissue)
# ������ �������� ����� �� ������
# counts[sigt,1]

# ������ ������ �������� �� �������� 
# p-val ������ ���� < 0.05 
sign_gene_age = filt_age$table$PValue < 0.05

# �� ���������� 312
sum(sign_gene_age) 

# ������ ��������  ����� �� ������
# counts[siga,1]

######------------------------------------------
#

#                    ����� 2

# �������������� ���� �������� ���� �� �� ������ �������
# ��� ������ ������������� ������������� (edgeR)���������� 1 �
# ����������� ���������� ��������) � 6 ���������.

#
######------------------------------------------

# �������� ����
sign_gene = sign_gene_age | sign_gene_tissue

# ����� �� ���������� 353, 
sum(sign_gene)

# ��� ������, ��� �� �����
# ����� �����, ������� ������� �� ����� ����������: 32

sign_norm_counts = t(norm_counts[sign_gene, ])

# ��������� ������� ����������
correlation_matrix = cor(sign_norm_counts, method = 'spearman')

# ������ �������� ������������� �������������
hcl_data = hclust(as.dist(1 - correlation_matrix))

# � ���������, ��� ����������
plot(hcl_data, hang = -1,  xlab = '����')
# � ����� �������, �� ������ ����� �� �����, ��-�� �������� ����������� ����� � ���������� ����������

# �������� �� 6 ���������
cl = cutree(hcl_data, k = 6)
rect.hclust(tree = hcl_data, k = 6, border = 1:6, cluster = cl)


######------------------------------------------
#

#                    ����� 3

# ������������ ���������� ������� ���� � �������� ���� �
# ��������� ���� (edgeR)z-score). ��������� ��� ������� ��������
# ����������� �������� z-score �� �������� ��� ����� ������

#
######------------------------------------------


z_score = scale(t(sign_norm_counts))
# �������� ��������
cerebellum = (tissue == 'C')
cortex = (tissue == "B")

dev.off() 
par(mfrow = c(3, 2))

#������ ��� ������� �������� ����������� �������� z-score �� �������� ��� ����� ������
for (i in 1:6) {
    a = colMeans(z_score[cl == i, ])
    plot(x = age[cerebellum], y = a[cerebellum], type='l', col = "#FF9466", lwd = 3,
         main = paste('�������', i), xlab = '�������, ����', ylab = '������� Z-score',
         ylim = c(min(a), max(a)))
    lines(x = age[cortex], y = a[cortex], type = 'l', col = '#75C5DC', lwd = 3)
}
legend(x=-17, y=4,
       legend=c("��������", "����"), fill=c("#FF9466", "#75C5DC"), xpd=NA, bty="n")

