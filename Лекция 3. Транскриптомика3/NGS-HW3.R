####################################
#
#  ДОМАШНЕЕ ЗАДАНИЕ 3. RNA-Seq
#
#       Горохов Никита
#
####################################


library(edgeR)

 
######------------------------------------------
#

#                    ЧАСТЬ 1

# Используя данные из предыдущего ДЗ при помощи edgeR)
# найти гены с межтканевыми и/или возрастными
# изменениями экспресии корректированное p-value < 0.05,
# межтканевые отличия должны быть не менее чем в два
# раза), возрастные изменения можно считать линейными по
# возрасту. Используйте модель ~ tissue + age

#
######------------------------------------------


# В count.out не отображаются названия экспериментов
# поэтому добавим их вручную

col_names = c("B14.5", "B15.5", "B17.5", "B20", "B34",
              "C14.5", "C15.5", "C17.5", "C20", "C34")

counts = read.table(file = 'count.out', stringsAsFactors = FALSE,
                    col.names = c("gene",col_names))

# Удалим гены, которые не экспрессируют во всех образцах
counts = counts[rowSums(counts[, 2:11]) >= 1, ] 

# В названиях присутстувет информация о типе ткани и возрасте
tissue = as.factor(substring(col_names, 1, 1))
age = as.numeric(substring(col_names, 2))

# Перейдем к обработке с помощью пакета edgeR
# сделаем по аналогии с лекциоными слайдами
edger = DGEList(counts = as.matrix(counts[,2:11]), group = tissue)
edger = calcNormFactors(edger, method = 'RLE') 
norm_counts = cpm(edger) 

# Будем использовать модель ~tissue + age
design = model.matrix(~ tissue + age)
edger = estimateGLMCommonDisp(edger,design)
edger = estimateGLMTrendedDisp(edger,design)
edger = estimateGLMTagwiseDisp(edger,design)
strict.disp = pmax(edger$tagwise.dispersion, edger$trended.dispersion, edger$common.dispersion)
plotBCV(edger)

#Теперь перейдем к построению и обучению GLM модели 
glm = glmFit(edger, design, dispersion=strict.disp)
lrt_tissue = glmLRT(glm, 2)
lrt_age = glmLRT(glm, 3)

# Отфильтруем данные с помощью topTags
filt_age = topTags(lrt_age, n = Inf, adjust.method = 'BH', sort.by = 'none')
filt_tissue = topTags(lrt_tissue, n = Inf, adjust.method = 'BH', sort.by = 'none')

# Теперь найдем значимые по типу ткани гены 
# p-val должно быть < 0.05 и межтканевые отличия 
# должны быть не менее чем в два раза
sign_gene_tissue = (filt_tissue$table$logFC >= 1 | filt_tissue$table$logFC <= -1) & (filt_tissue$table$PValue < 0.05)

#Их количество: 73
sum(sign_gene_tissue)
# Список значимых генов по тканям
# counts[sigt,1]

# Теперь найдем значимые по возрасту 
# p-val должно быть < 0.05 
sign_gene_age = filt_age$table$PValue < 0.05

# Их получилось 312
sum(sign_gene_age) 

# Списко значимых  генов по тканям
# counts[siga,1]

######------------------------------------------
#

#                    ЧАСТЬ 2

# Скалстеризуйте гены значимые хотя бы по одному фактору
# при помощи иерархической кластеризации (edgeR)расстояние 1 —
# коэффициент корреляции Спирмана) в 6 кластеров.

#
######------------------------------------------

# Обединим гены
sign_gene = sign_gene_age | sign_gene_tissue

# Итого их получилось 353, 
sum(sign_gene)

# что меньше, чем их сумма
# Всего генов, которые значимы по обоим категориям: 32

sign_norm_counts = t(norm_counts[sign_gene, ])

# Посчитаем матрицу корреляции
correlation_matrix = cor(sign_norm_counts, method = 'spearman')

# Теперь применим иерархическую кластеризацию
hcl_data = hclust(as.dist(1 - correlation_matrix))

# И посмотрим, что получилось
plot(hcl_data, hang = -1,  xlab = 'Гены')
# В целом понятно, но осбого смсла не несет, из-за большого количестова генов и маленького разрешения

# Разобьем на 6 кластеров
cl = cutree(hcl_data, k = 6)
rect.hclust(tree = hcl_data, k = 6, border = 1:6, cluster = cl)


######------------------------------------------
#

#                    ЧАСТЬ 3

# Отшкалируйте экспрессию каждого гена к среднему ноль и
# дисперсии один (edgeR)z-score). Нарисуйте для каждого кластера
# зависимость среднего z-score от возраста для обоих тканей

#
######------------------------------------------


z_score = scale(t(sign_norm_counts))
# Выбираем мозжечок
cerebellum = (tissue == 'C')
cortex = (tissue == "B")

dev.off() 
par(mfrow = c(3, 2))

#Рисуем для каждого кластера зависимость среднего z-score от возраста для обоих тканей
for (i in 1:6) {
    a = colMeans(z_score[cl == i, ])
    plot(x = age[cerebellum], y = a[cerebellum], type='l', col = "#FF9466", lwd = 3,
         main = paste('Кластер', i), xlab = 'Возраст, день', ylab = 'Средний Z-score',
         ylim = c(min(a), max(a)))
    lines(x = age[cortex], y = a[cortex], type = 'l', col = '#75C5DC', lwd = 3)
}
legend(x=-17, y=4,
       legend=c("Мозжечок", "Кора"), fill=c("#FF9466", "#75C5DC"), xpd=NA, bty="n")

