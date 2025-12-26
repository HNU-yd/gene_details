# Ensembl/GFF/VCF/FASTA Genome Stats Tree — 输出字典字段说明

本项目会将 `data/<species>/` 下的 FASTA（参考序列）、GFF/GFF3（注释）与 VCF（变异）整合，输出分层树状统计结构，并支持：
- 每条染色体单独输出（人类常见：GFF/VCF 按染色体拆分）
- 最终输出一个物种级汇总统计（不包含所有 gene 子树，避免巨大文件）

---

## 1. 输出目录与文件

输出根目录固定为：

```

./stats/<species>/

````

包含文件：

- `chrom_<chrom>.tree.json.gz`  
  单条染色体的完整树结构（包含 gene → exon → cds/utr/intron 等细分段）与统计字段。

- `chrom_<chrom>.nodes.tsv.gz`  
  染色体树的扁平化节点表（便于 pandas / grep 分析）。

- `species.total.summary.json.gz`  
  物种级汇总（由多个染色体的“浅拷贝 summary 节点”组成，不包含 gene 级子树）。

- `species.total.summary.nodes.tsv.gz`  
  物种级汇总的扁平化节点表。

---

## 2. JSON 结构总览

每个 `*.tree.json.gz` 解压后是一个 JSON 对象，对应一个根节点（Node）：

- 染色体文件：根节点 `type="chromosome"`
- 物种汇总文件：根节点 `type="species"`

Node 结构如下：

```json
{
  "type": "...",
  "id": "...",
  "chrom": "...",
  "start": ...,
  "end": ...,
  "length_span": ...,
  "coverage": [[s1,e1],[s2,e2],...],
  "covered_length": ...,
  "variant_count": ...,
  "missing_N_count": ...,
  "missing_N_ratio": ...,
  "missing_N_runs": ...,
  "child_quant": {...},
  "desc_leaf_quant": {...},
  "children": [...]
}
````

---

## 3. Node 字典字段逐项解释

### 3.1 标识类字段

* `type`（string）
  节点类型，决定它在树中的层级角色。常见类型见「第 4 节」。

* `id`（string）
  节点唯一标识符（同一输出文件内唯一）。
  常见示例：

  * `chrom:1`
  * `chrom:1:gene_union_group`
  * `gene:Os01g0100200`
  * `gene:Os01g0100200:exon_group`
  * `...:cds_group:seg:3`

* `chrom`（string | null）
  规范化后的染色体名（如 `"1"`, `"2"`, `"X"`）。
  物种根节点 `species` 通常为 `null`。

---

### 3.2 坐标与区间字段（1-based 闭区间）

> 坐标体系：**1-based，闭区间 [start, end]**
> 与 GFF/VCF 常规坐标一致。

* `start`（int | null）
  节点覆盖区域的起点（闭区间）。

* `end`（int | null）
  节点覆盖区域的终点（闭区间）。

* `length_span`（int）
  节点跨度长度（包围盒跨度）：
  `length_span = end - start + 1`
  不考虑中间空洞。

* `coverage`（array of [int,int]）
  节点的**实际覆盖区间集合**，为若干不重叠的闭区间，通常已经合并（merged）。
  例：exon_union 可为多段：

  ```json
  "coverage": [[11218,12060],[12152,12435]]
  ```

* `covered_length`（int）
  coverage 内所有区间长度之和：
  `covered_length = Σ (ei - si + 1)`
  对多段结构（exon_union、gene_union 等）比 `length_span` 更真实。

---

### 3.3 变异统计字段（来自 VCF）

* `variant_count`（int）
  节点所覆盖区域内的变异数量。

#### 重要：`variant_count` 有两类语义

**A. Unique 计数（不等于子节点之和）**
以下节点的 `variant_count` 是“按 VCF 记录唯一计数”，不会强制等于子节点累加：

* `chromosome`：该染色体 VCF 中记录数（只计本 chrom 的行）
* `gene_union_group`：落在“基因并集区间”内的 VCF 记录数（唯一计数）
* `non_gene_group`：落在“非基因区间”内的 VCF 记录数（唯一计数）

原因：一个变异可能与多个 gene 或多个 segment 重叠；若做子节点求和会重复计数。

**B. Sum-of-leaves 计数（按叶子累加）**
除上述 unique 节点外，多数中间节点（如 `gene`、`exon_group`、`cds_group` 等）会将
`variant_count` 聚合为“其子树叶子节点 `variant_count` 的总和”。

这意味着：

* `gene.variant_count` 可能 **大于** 该 gene 区间内“唯一变异数”
* 适用于分析“结构分布”（变异落在 CDS/UTR/INTRON 的相对量），但不等价于唯一位点数

---

### 3.4 缺失（N）统计字段（来自 FASTA）

缺失定义：FASTA 序列中的字符 `N`（连续 N-run 也被记录）。
统计范围：默认在该节点的 **coverage** 覆盖范围内计算（不是 length_span）。

* `missing_N_count`（int）
  coverage 覆盖范围内 `N` 碱基数量总和。

* `missing_N_ratio`（float）
  `N` 占 coverage 的比例：
  `missing_N_ratio = missing_N_count / covered_length`
  若 `covered_length=0` 则为 0。

* `missing_N_runs`（int）
  coverage 覆盖到的 N-run 命中次数（可理解为“与多少段连续 N 区间发生交叠”）。

---

### 3.5 下一级量化统计 `child_quant`

`child_quant` 用于满足需求：
“每一级必须量化下一级：数量、最小长度、最大长度、均值、方差、缺失比例统计、变异总数”。

结构：

```json
"child_quant": {
  "by_type": {
    "<child_type>": {
      "count": ...,
      "len_min": ...,
      "len_max": ...,
      "len_mean": ...,
      "len_var": ...,
      "missing_ratio_min": ...,
      "missing_ratio_max": ...,
      "missing_ratio_mean": ...,
      "missing_ratio_var": ...,
      "variant_sum": ...
    }
  }
}
```

字段解释：

* `by_type`：按“直接子节点类型”分组。

每个 `<child_type>` 下：

* `count`：该类型**直接子节点**数量
  例：`exon_segment` 的子节点通常有 `cds_group / utr_group / exon_other_group` 各一个。

* `len_min / len_max / len_mean / len_var`：
  以子节点的 `covered_length` 为样本计算统计量

  * `len_var` 是总体方差（population variance）

* `missing_ratio_*`：
  以子节点的 `missing_N_ratio` 为样本计算统计量

* `variant_sum`：
  该类型所有子节点 `variant_count` 的直接求和
  注意：若子节点是 unique 类型，则求和是 unique 的；否则是结构累加意义。

---

### 3.6 子树叶子汇总 `desc_leaf_quant`

`desc_leaf_quant` 用于在任意层级直接回答：

* “该节点下面 CDS/UTR/INTRON/NON_GENE 等叶子分别有多少段？”
* “这些叶子段分别累计了多少变异（结构累加）？”

结构：

```json
"desc_leaf_quant": {
  "by_leaf_type": {
    "cds_segment": {"count": ..., "variant_sum": ...},
    "utr_segment": {"count": ..., "variant_sum": ...},
    "intron_segment": {"count": ..., "variant_sum": ...},
    "exon_other_segment": {"count": ..., "variant_sum": ...},
    "non_gene_segment": {"count": ..., "variant_sum": ...}
  }
}
```

字段解释：

* `count`：该叶子类型在该节点子树中的叶子段数量
* `variant_sum`：这些叶子段的 `variant_count` 总和（结构分布意义，不保证唯一）

---

### 3.7 子节点列表 `children`

* `children`（array of Node）
  节点的直接子节点列表。树形组织见下一节。

---

## 4. 树结构与节点类型（type）

### 4.1 染色体树结构（chrom_<chrom>.tree.json.gz）

```
chromosome
├── gene_union_group
│   └── gene (many)
│       ├── exon_group
│       │   └── exon_segment (many)
│       │       ├── cds_group -> cds_segment*        (many)
│       │       ├── utr_group -> utr_segment*        (many)
│       │       └── exon_other_group -> exon_other_segment* (many)
│       └── non_exon_group -> intron_segment*        (many)
└── non_gene_group -> non_gene_segment*              (many)
```

### 4.2 type 语义说明

* `chromosome`
  染色体整体节点，coverage 通常是 `[1, chrom_len]`。

* `gene_union_group`
  染色体上所有 gene span 的并集（多个 `[gene.start, gene.end]` 合并得到）。

* `gene`
  单个基因节点，coverage 通常为 gene span `[start,end]`。

* `exon_group`
  该 gene 内所有 transcript 的 exon 合并并集（exon_union）。

* `exon_segment`
  exon_union 的一个连续区段。

* `cds_group` / `utr_group` / `exon_other_group`
  对单个 `exon_segment` 的进一步拆分：

  * `cds_group`：CDS union 与该 exon_segment 的交集
  * `utr_group`：UTR union 与该 exon_segment 的交集
  * `exon_other_group`：不属于 CDS/UTR 的剩余部分

* `cds_segment` / `utr_segment` / `exon_other_segment`
  叶子段（用于结构统计与变异分布累计）。

* `non_exon_group`
  gene span 中减去 exon_union 的部分（intron 集合）。

* `intron_segment`
  intron 的叶子段。

* `non_gene_group`
  染色体中不属于 gene_union 的区域（gene_union 的补集）。

* `non_gene_segment`
  非基因区的叶子段。

---

## 5. TSV 输出字段说明（chrom_<chrom>.nodes.tsv.gz 等）

每行一个节点（前序遍历），字段：

* `depth`：节点深度（根为 0）
* `parent_id`：父节点 id（根为空）
* `id`：节点 id
* `type`：节点类型
* `chrom`：染色体
* `start`, `end`：坐标（1-based inclusive）
* `length_span`：跨度长度（包围盒）
* `covered_length`：coverage 长度和
* `variant_count`：变异计数（注意 unique vs sum-of-leaves）
* `missing_N_count`：N 总数（coverage 范围）
* `missing_N_ratio`：N 占比（coverage 范围）
* `missing_N_runs`：命中 N-run 次数（coverage 范围）

---

## 6. 常用查询方式（面向你的需求）

### 6.1 “exon 下面分别有多少 CDS / UTR？”

看 `exon_segment` 节点：

* 直接子类型数量（group 级）：

  * `child_quant.by_type.cds_group.count`
  * `child_quant.by_type.utr_group.count`

如果你要叶子段数量（segment 级）：

* `desc_leaf_quant.by_leaf_type.cds_segment.count`
* `desc_leaf_quant.by_leaf_type.utr_segment.count`

### 6.2 “CDS 下有多少变异？UTR 下有多少变异？”

看任意层级节点：

* `desc_leaf_quant.by_leaf_type.cds_segment.variant_sum`
* `desc_leaf_quant.by_leaf_type.utr_segment.variant_sum`

说明：这是“按叶子段累加”的结构分布计数，并不保证唯一位点去重。

---

## 7. 人类按染色体拆分文件的适配逻辑

当输入是“每条染色体一个 GFF/VCF 文件”时，程序会：

* 对每个 GFF：

  * 优先读取 `##sequence-region <seqid> ...` 推断 chrom
  * 否则使用第一条 feature 行的 seqid（第 1 列）推断 chrom

* 对每个 VCF：

  * 使用第一条非 header 的变异行（第 1 列 `CHROM`）推断 chrom

然后逐条染色体处理并输出 `chrom_<chrom>.*`，最后写 `species.total.summary.*`。

---

```
```
