# Monocle 轨迹生成算法细节 (基于 DDRTree)

本文档详细阐述了 Monocle 2.x 版本中，在 DDRTree 降维之后如何生成细胞轨迹的算法细节。

**关于 UMAP 降维的说明：**

Monocle 2.x 版本（如本分析所基于的代码库）的 `reduceDimension` 函数默认支持的降维方法包括 DDRTree、ICA、tSNE 等，但**不直接内置 UMAP 降维功能**。如果您在 Monocle 外部使用 UMAP 进行了降维，并将降维后的坐标导入 Monocle，那么 Monocle 仍然可以基于这些坐标构建轨迹。然而，Monocle 内部的轨迹学习算法（特别是伪时间计算和状态分配）是针对其内置的降维方法（尤其是 DDRTree）进行优化的。

以下将详细介绍 Monocle 中 DDRTree 降维和轨迹生成的算法流程。

## 1. 降维 (Dimensionality Reduction) - `reduceDimension` 函数

当 `reduceDimension` 函数的 `reduction_method` 参数设置为 `"DDRTree"` 时，Monocle 会使用 Discriminative Dimensionality Reduction with Trees (DDRTree) 算法进行降维。

**DDRTree 算法概述：**

DDRTree 是一种流形学习算法，旨在发现高维数据中的低维非线性结构。它通过学习一个“主图”（principal graph）来表示数据的拓扑结构。这个主图由一系列节点（代表细胞状态的中心点）和连接这些节点的边组成，这些节点和边共同构成了细胞在基因表达空间中演变轨迹的骨架。

**算法步骤：**

1.  **数据准备：**
    *   输入：经过归一化（例如 `log` 或 `vstExprs`）和可选的批次效应移除后的基因表达矩阵 `FM`。
    *   对 `FM` 进行转置和标准化（如果 `scaling` 为 `TRUE`），确保数据适合 DDRTree 算法的输入要求。

2.  **调用 DDRTree 核心：**
    *   Monocle 内部调用 `DDRTree` R 包中的 `DDRTree` 函数。
    *   该函数接收处理后的表达矩阵 `FM` 和目标降维维度 `max_components` 作为输入。
    *   如果 `auto_param_selection` 为 `TRUE` 且细胞数量较大，Monocle 会自动计算 `ncenter` 参数（主图中的节点数量），以优化 DDRTree 的性能。

3.  **DDRTree 输出：**
    *   `ddrtree_res$W`：权重矩阵，用于将高维数据映射到低维空间。存储在 `cds@reducedDimW`。
    *   `ddrtree_res$Z`：降维后的细胞坐标（latent positions of cells），表示每个细胞在低维空间中的位置。存储在 `cds@reducedDimS`。
    *   `ddrtree_res$Y`：主图的节点坐标（principal graph nodes），表示轨迹骨架上的关键点。存储在 `cds@reducedDimK`。

4.  **构建最小生成树 (MST)：**
    *   基于 `ddrtree_res$Y`（主图节点坐标）计算节点间的欧氏距离。
    *   使用 `igraph` 包的 `graph.adjacency` 和 `minimum.spanning.tree` 函数，构建连接这些主图节点的最小生成树 (MST)。这个 MST 代表了学习到的细胞轨迹的拓扑结构。存储在 `cds@minSpanningTree`。

5.  **细胞投影到 MST：**
    *   `findNearestPointOnMST` 函数：将每个细胞（在 `reducedDimS` 中）投影到 MST 上最近的节点。
    *   `project2MST` 函数：进一步将每个细胞投影到 MST 上最近的线段上，得到更精确的投影点 `P`。这些投影点将用于后续的伪时间计算。
    *   更新 `cds@cellPairwiseDistances` 和 `cds@minSpanningTree` 以反映这些投影。

## 2. 轨迹排序 (Trajectory Ordering) - `orderCells` 函数

`orderCells` 函数在降维的基础上，为每个细胞分配伪时间 (Pseudotime) 和状态 (State)，从而构建完整的细胞轨迹。

**算法步骤：**

1.  **选择根细胞 (Root Cell)：**
    *   `select_root_cell` 函数负责确定轨迹的起始点。
    *   如果用户未指定 `root_state`，Monocle 会默认选择 MST 直径（最长路径）的第一个端点作为根细胞。
    *   如果指定了 `root_state`，则从属于该状态的细胞中选择一个合适的细胞作为根。
    *   对于 DDRTree 降维，根细胞会被进一步映射到主图上最近的节点。

2.  **基于 DDRTree 的轨迹提取 (`extract_ddrtree_ordering`)：**
    *   **MST 遍历：** 从选定的根细胞开始，对 `cds@minSpanningTree`（连接主图节点的 MST）进行深度优先搜索 (DFS) 遍历。
    *   **伪时间计算：**
        *   根细胞的伪时间被设置为 0。
        *   对于 MST 上的每个细胞，其伪时间是其到根细胞在 MST 上的最短路径距离。这个距离是基于细胞在降维空间中的欧氏距离累积得到的。
        *   `dp[curr_node_name, parent_node_name]` 表示当前细胞与其父细胞在降维空间中的距离，这个距离被累加到父细胞的伪时间上，得到当前细胞的伪时间。
    *   **状态分配：**
        *   Monocle 根据 MST 上的分支点（即度数大于 2 的节点）来识别不同的细胞状态或分支。
        *   当 DFS 遍历遇到一个分支点时，会创建一个新的细胞状态。
        *   所有从该分支点开始沿特定路径延伸的细胞，都会被分配到同一个状态，直到遇到下一个分支点或轨迹终点。
        *   `curr_state` 变量用于跟踪当前的状态编号，并在遇到分支时递增。

3.  **更新 CellDataSet 对象：**
    *   计算出的伪时间 (`pseudo_time`) 和细胞状态 (`cell_state`) 会被添加到 `cds` 对象的 `pData`（表型数据）中，作为 `Pseudotime` 和 `State` 列。
    *   `cds@auxOrderingData` 中会存储轨迹相关的辅助信息，例如根细胞和分支点。

## 总结

Monocle 通过 DDRTree 算法将高维单细胞基因表达数据降维到低维空间，并学习一个代表细胞状态连续变化的主图。随后，通过在主图上构建最小生成树，并从一个选定的根细胞开始遍历该树，Monocle 能够计算每个细胞的伪时间，并根据轨迹中的分支点识别不同的细胞状态。这一系列步骤共同揭示了细胞分化或发育的动态过程。
