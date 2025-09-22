1. 功能富集分析（最核心的分析）
GO富集分析
你的基因名是形如 10123_at、2487_at 的探针 ID（Affymetrix 芯片格式），
而不是标准的 HGNC 基因符号（如 TP53、IL6），所以 mygene 无法识别它们。这种情况非常常见于 GEO 芯片数据分析
一般要下载GPLannota解释文件，然后转换，再映射会官方基因名
例如
✅ 解决方案：先将探针 ID 映射为基因符号
你需要一个 GEO 平台注释文件（如 GPL570），然后在 Python 中完成映射。下面是完整流程：

🧬 Step 1：读取注释文件并构建映射
```python
import pandas as pd

# 读取 GPL570 注释文件（你需提前下载并解压）
annot = pd.read_csv('GPL570.annot', sep='\t', skiprows=16)
annot = annot[['ID', 'Gene Symbol']].dropna()

# 构建 probe → symbol 映射字典
probe_to_symbol = dict(zip(annot['ID'], annot['Gene Symbol']))

# 应用到 deg_df
deg_df['Gene_Symbol'] = deg_df.index.map(probe_to_symbol)

# 清洗多重注释（如 "TP53 /// MDM2"）
deg_df['Gene_Symbol'] = deg_df['Gene_Symbol'].str.split(' /// ').str[0]
```
✅ Step 2：重新运行你的分析代码
```python
# DEG 上调基因（log2FC < -1）
deg_up = deg_df[(deg_df['significant']) & (deg_df['log2fc'] < -1)].sort_values('fdr').head(10)
gene_symbols = deg_up['Gene_Symbol'].dropna().unique().tolist()
```

---
KEGG通路分析
```python
# KEGG通路富集
kegg_results = gp.enrichr(gene_list=deg_genes,
                         gene_sets=['KEGG_2021_Human'],
                         organism='Human')

# 通路可视化
plt.figure(figsize=(12, 8))
kegg_results.results.head(15).plot.bar(x='Term', y='Adjusted P-value')
plt.title('KEGG Pathway Enrichment')
plt.xticks(rotation=45, ha='right')
```

---
2. 蛋白互作网络（PPI）分析
使用STRING数据库
```python
import networkx as nx
import requests
#step1 构建PPI网络
#step2 构建网络
#step3 识别Hub基因
```

---
3. 转录因子调控分析
```python
# 使用ChEA3或TRRUST数据库分析转录因子
def tf_enrichment_analysis(deg_genes):
    # 这里可以使用gseapy的转录因子数据库
    tf_results = gp.enrichr(gene_list=deg_genes,
                           gene_sets=['ChEA_2022', 'TRRUST_Transcription_Factors_2019'],
                           organism='Human')
    return tf_results

# 识别关键的转录因子
tf_results = tf_enrichment_analysis(deg_genes)
```

---
4. 药物-基因相互作用分析
```python
# 药物重定位分析
def drug_enrichment_analysis(deg_genes):
    drug_results = gp.enrichr(gene_list=deg_genes,
                             gene_sets=['DrugMatrix', 'Drug_Perturbations_from_GEO_2022'],
                             organism='Human')
    return drug_results

# 寻找潜在的治疗药物
drug_results = drug_enrichment_analysis(deg_genes)
```

---
5. 机器学习与生物标志物识别
```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

# 使用差异基因构建预测模型
def build_prediction_model(expression_data, deg_genes, labels):
    # 选择差异基因的表达数据
    X = expression_data.loc[deg_genes].T
    y = labels
    
    # 训练随机森林模型
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # 评估模型
    y_pred = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_pred)
    
    # 获取基因重要性
    feature_importance = pd.DataFrame({
        'gene': deg_genes,
        'importance': model.feature_importances_
    }).sort_values('importance', ascending=False)
    
    return model, auc, feature_importance
```

---
6. 生存分析（如果有临床数据）

---
7.与其他组学数据库结合
```python
# 与甲基化数据整合
def integrate_methylation(deg_genes, methylation_data):
    # 寻找差异基因的甲基化状态
    integrated_results = []
    
    for gene in deg_genes:
        if gene in methylation_data.index:
            gene_methylation = methylation_data.loc[gene]
            # 分析表达与甲基化的相关性
            # ...
    
    return integrated_results

# 与蛋白组数据整合
def integrate_proteomics(deg_genes, proteomics_data):
    # 检查mRNA-protein表达一致性
    concordance_analysis = {}
    
    for gene in deg_genes:
        if gene in proteomics_data.index:
            mrna_expression = expression_data.loc[gene]
            protein_expression = proteomics_data.loc[gene]
            correlation = np.corrcoef(mrna_expression, protein_expression)[0, 1]
            concordance_analysis[gene] = correlation
    
    return concordance_analysis
```

---
8. 可视化与结果展示
```text
 # 1. 火山图
# 2. GO富集气泡图
 # 3. PPI网络图
 # 4. 通路图
    # 可以添加KEGG通路图
#5. 热图
#6.PCA图
总的来说最好每个分析都能出张图
```
