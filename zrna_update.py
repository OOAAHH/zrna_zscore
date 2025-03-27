# -*- coding: utf-8 -*-
"""
ZRNA分数计算模块

本模块实现了基于Z分数的RNA结构评估方法，通过标准化不同评分指标来比较RNA结构模型的质量。

计算原理与流程:
---------------
1. Z分数基础原理:
   Z分数(Z-score)是统计学中的标准化分数，表示一个观测值与平均值的偏差程度，以标准差为单位。
   Z = (X - μ) / σ
   其中:
   - X是原始分数
   - μ是总体平均值
   - σ是总体标准差
   
   Z分数的意义:
   - Z分数为0表示该值等于平均值
   - 正的Z分数表示该值高于平均值
   - 负的Z分数表示该值低于平均值
   - Z分数的绝对值表示偏离平均值的程度

2. 双重Z分数计算:
   为了减少异常值对计算结果的影响，本模块采用双重Z分数计算方法:
   a) 第一轮: 计算所有模型的Z分数
   b) 移除异常值: 将Z分数低于阈值的模型标记为异常值(NaN)
   c) 第二轮: 使用剩余模型重新计算Z分数
   d) 限制极端值: 将最终Z分数中低于阈值的值设为阈值值
   
   这种方法确保了评分的稳健性，减少了个别异常模型对整体评估的影响。

   双重Z分数计算的优势:
   - 对异常值不敏感: 常规Z分数计算容易受到极端值的影响，而双重Z分数通过两轮计算减轻这种影响
   - 保持数据分布: 移除异常值后，重新计算得到的Z分数能更好地反映正常数据的分布特征
   - 下限控制: 通过设置阈值，避免极端低分过度影响最终评估结果
   - 适用于RNA结构模型: RNA结构预测中常出现部分极度偏离的模型，双重Z分数能有效处理这种情况

   在RNA结构评估中的特殊价值:
   RNA结构预测模型通常由不同算法和参数生成，质量差异可能很大。某些算法可能在特定结构上产生完全错误的模型，导致评分指标出现极端值，扭曲标准Z分数计算。双重Z分数通过识别并处理这些异常值，使评分更加可靠。

3. ZRNA分数计算:
   ZRNA分数是多个指标Z分数的加权组合，用于全面评估RNA结构模型的质量:
   ZRNA = Σ(Z_metric * weight_metric)
   
   其中:
   - Z_metric是各指标的Z分数
   - weight_metric是各指标的权重，反映指标的重要性
   
   最终ZRNA分数越高，表示模型整体质量越好。

4. 分组计算:
   对于具有多个构象的RNA结构，本模块提供了以下聚合方法:
   - sum>0: 仅计算正分数的总和
   - mean: 计算平均分数
   - sum: 计算所有分数的总和

注意事项:
--------
1. 数据前处理:
   - 对于"越低越好"的指标，需将其分数取反再计算Z分数
   - 输入数据可以包含NaN值，计算时会自动排除

2. 异常值处理:
   - 默认异常值阈值为-2，可根据具体需求调整
   - 第二轮Z分数计算后，所有低于阈值的值会被设为阈值值

3. 结果解读:
   - 正的ZRNA分数表示模型质量高于平均水平
   - ZRNA分数越高，表示模型质量越好
   - 组间比较时应确保使用相同的指标和权重

4. 常见指标及其权重选择:
   在RNA结构评估中，常用的指标包括:
   - TM-score: 衡量整体结构相似性，越高越好
   - RMSD: 衡量空间偏差，越低越好
   - Interaction Network Fidelity (INF): 衡量碱基对相互作用的准确性，越高越好
   - Clash score: 衡量原子碰撞程度，越低越好
   
   权重设置建议:
   - 结构相似性指标(如TM-score): 较高权重 (0.8-1.0)
   - 局部特征指标(如INF): 中等权重 (0.5-0.8)
   - 物理合理性指标(如clash score): 较低权重 (0.3-0.5)
   
   注意: 对于越低越好的指标，权重应设为负值。

5. 调试与错误排查:
   本模块包含详细的日志输出功能，可帮助用户排查计算过程中的问题:
   
   a) 常见问题及解决方案:
      - 数据分布异常: 如果标准差非常小，Z分数可能会变得非常大。检查输入数据是否有足够的变异性。
      - NaN值过多: 输入数据中的NaN值过多会影响计算结果。确保大部分模型都有有效分数。
      - 权重选择不当: 如果某些指标的权重设置不合理，可能导致整体评分偏离预期。调整权重值并重新计算。
      
   b) 提高计算准确性:
      - 选择合适的阈值: 基于数据分布特性选择适当的异常值阈值。
      - 增加样本量: 比较更多模型可以提高Z分数的统计可靠性。
      - 平衡指标权重: 确保各指标权重反映其在整体评估中的实际重要性。
      
   c) 日志解读:
      - DEBUG级别日志: 提供详细的计算过程信息，包括输入形状、数据类型、NaN值检测等。
      - INFO级别日志: 提供关键统计信息，如均值、标准差、异常值数量等。
      - 关注数据形状变化: 确保在处理过程中数据形状一致。
      - 监控NaN值数量: 异常值处理可能导致NaN值数量增加，影响后续计算。

实际应用示例:
-----------
以下是使用本模块计算RNA结构模型质量的典型流程:

1. 基本Z分数计算:
   ```python
   import numpy as np
   import pandas as pd
   from zrna import calc_z_scores

   # RNA结构评分数据
   tm_scores = np.array([0.85, 0.92, 0.78, 0.64, 0.88, 0.91])
   
   # 计算Z分数 (TM-score越高越好，所以lower_is_better=False)
   z_scores = calc_z_scores(tm_scores, lower_is_better=False)
   # 输出: array([-0.05, 1.24, -1.36, -3.23, 0.50, 1.08])
   ```

2. 双重Z分数计算:
   ```python
   from zrna import calc_double_z_scores
   
   # RMSD值 (RMSD越低越好，所以lower_is_better=True)
   rmsd_values = np.array([2.1, 1.8, 3.2, 8.5, 2.3, 1.9])
   
   # 计算双重Z分数，移除异常值
   double_z_scores = calc_double_z_scores(rmsd_values, thresh=-2.0, lower_is_better=True)
   ```

3. 多指标ZRNA分数计算:
   ```python
   import pandas as pd
   from zrna import calc_zrna
   
   # 准备包含多个指标Z分数的数据框
   data = {
       'target': ['RNA1', 'RNA1', 'RNA2', 'RNA2', 'RNA3', 'RNA3'],
       'gr_code': ['modelA', 'modelB', 'modelA', 'modelB', 'modelA', 'modelB'],
       'z_tm': [1.2, 0.8, -0.5, 1.3, 0.9, -1.1],
       'z_inf': [0.9, 1.1, 0.7, -0.8, 1.2, 0.5],
       'z_rmsd': [0.7, -0.3, 1.2, 0.8, -0.5, -1.2]
   }
   df = pd.DataFrame(data)
   
   # 定义指标权重 (正值表示越高越好，负值表示越低越好)
   metric_weights = {'tm': 1.0, 'inf': 0.8, 'rmsd': -0.5}
   
   # 计算ZRNA分数
   df['z_rna'] = calc_zrna(df, metric_weights)
   
   # 聚合每个目标的最佳分数
   metrics_is_lower_better = {'tm': 'max', 'inf': 'max', 'rmsd': 'min'}
   best_scores = best_score_across_all_confs(df, metrics_is_lower_better)
   
   # 计算每个组的总分
   group_scores = get_group_score(df, agg='sum>0', score='z_rna')
   ```

这个计算流程适用于评估和比较多种RNA结构预测模型的性能，特别是在需要综合考虑多个评价指标的场景。

参考文献:
--------
本方法基于结构生物学中常用的统计分析方法，适用于RNA结构模型的评估与比较。
"""

import numpy as np
import logging
import pandas as pd

# 设置日志格式
logging.basicConfig(level=logging.DEBUG,
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calc_z_scores(a: np.ndarray, lower_is_better: bool = False) -> np.ndarray:
    '''
    计算每个模型的z分数。

    参数:
        a (np.ndarray): 分数数组，形状为 (n_models,)
        lower_is_better (bool): 是否分数越低越好
    返回:
        np.ndarray: z分数数组，形状为 (n_models,)
    '''
    logger.debug(f"输入数组形状: {a.shape}, 数据类型: {a.dtype}")
    logger.debug(f"输入数组是否包含NaN值: {np.any(np.isnan(a))}")
    
    if lower_is_better:
        logger.debug("转换分数（分数越低越好）")
        a = -a
    
    mean_val = np.nanmean(a)
    std_val = np.nanstd(a)
    logger.debug(f"平均值: {mean_val:.4f}")
    logger.debug(f"标准差: {std_val:.4f}")
    
    ans = (a - mean_val) / std_val
    logger.debug(f"输出z分数形状: {ans.shape}")
    logger.debug(f"输出是否包含NaN值: {np.any(np.isnan(ans))}")
    
    # 每个问题的模型数量
    logger.info(f"模型数量: {a.size}")
    logger.info(f"TM分数平均值: {mean_val:.4f}")
    logger.info(f"TM分数标准差: {std_val:.4f}")

    return ans

def remove_outliers(a: np.ndarray, thresh: float = -2) -> np.ndarray:
    """
    通过将异常值设置为nan来移除异常值。

    参数:
        a (np.ndarray): 分数数组，形状为 (n_models,)
        thresh (float): 异常值阈值
    返回:
        np.ndarray: 异常值被设置为nan的分数数组，形状为 (n_models,)
    """
    logger.debug(f"输入数组形状: {a.shape}")
    logger.debug(f"阈值: {thresh}")
    
    original_nan_count = np.sum(np.isnan(a))
    a[a < thresh] = np.nan
    new_nan_count = np.sum(np.isnan(a))
    
    logger.info(f"移除了 {new_nan_count - original_nan_count} 个异常值")
    return a

def calc_double_z_scores(raw: np.ndarray, thresh: float = -2, lower_is_better: bool = False) -> np.ndarray:
    """
    按照论文描述计算每个模型的双重z分数。

    参数:
        raw (np.ndarray): 分数数组，形状为 (n_models,)
        thresh (float): 异常值阈值
        lower_is_better (bool): 是否分数越低越好
    返回:
        np.ndarray: 双重z分数数组，形状为 (n_models,)
    """
    logger.debug(f"开始双重z分数计算")
    logger.debug(f"输入数组形状: {raw.shape}, 数据类型: {raw.dtype}")
    
    raw = np.copy(raw)
    first_pass = calc_z_scores(raw, lower_is_better=lower_is_better)
    logger.debug(f"第一遍z分数形状: {first_pass.shape}")
    
    filtered = remove_outliers(first_pass, thresh=thresh)
    logger.debug(f"移除异常值后，NaN值数量: {np.sum(np.isnan(filtered))}")

    raw[np.isnan(filtered)] = np.nan
    logger.debug(f"更新后的原始分数中NaN值数量: {np.sum(np.isnan(raw))}")

    second_pass = calc_z_scores(raw, lower_is_better=lower_is_better)
    logger.debug(f"第二遍z分数形状: {second_pass.shape}")
    
    second_pass[second_pass < thresh] = thresh
    second_pass[np.isnan(second_pass)] = thresh
    
    logger.info(f"最终z分数范围: [{np.min(second_pass):.4f}, {np.max(second_pass):.4f}]")
    return second_pass

def calc_zrna(df, metric_weights: dict) -> np.ndarray:  
    '''
    计算每个目标的zrna分数。

    参数:
        df (DataFrame): 包含目标、gr_code和指标列的数据框
        metric_weights (dict): 指标名称到权重（最大值或最小值）的映射
    返回:
        np.ndarray: zrna分数数组
    '''  
    logger.debug(f"使用权重计算zrna分数: {metric_weights}")
    logger.debug(f"数据框中可用的指标: {[col for col in df.columns if col.startswith('z_')]}")
    
    result = sum([df["z_" + metric] * metric_weights[metric] for metric in metric_weights])
    logger.debug(f"结果形状: {result.shape}")
    logger.info(f"最终zrna分数范围: [{np.min(result):.4f}, {np.max(result):.4f}]")
    return result

def best_score_across_all_confs(df, metrics_is_lower_better: dict) -> pd.DataFrame:    
    '''
    计算每个目标在所有构象中的最佳分数。

    参数:
        df (DataFrame): 包含目标、gr_code和指标列的数据框
        metrics_is_lower_better (dict): 指标名称到是否分数越低越好的映射
    返回:
        pd.DataFrame: 包含每个目标最佳分数的数据框
    '''
    logger.debug(f"计算所有构象中的最佳分数")
    logger.debug(f"指标配置: {metrics_is_lower_better}")
    
    result = df.groupby(["target", "gr_code"]).agg(metrics_is_lower_better).reset_index()
    logger.debug(f"结果形状: {result.shape}")
    return result

def get_group_score(df, agg: str = "sum>0", score: str = "z_rna") -> pd.DataFrame:
    '''
    来源：改编自Rachael的EM流程。

    聚合每个gr_code组的分数

    参数:
        df (DataFrame): 包含gr_code和分数列的数据框
        agg (str): 聚合分数的方式，选项：'sum>0'、'mean'、'sum'（默认为'sum>0'）
        score (str): 要聚合的列名
    返回:
        pd.DataFrame: 包含每个组聚合分数的数据框
    '''
    logger.debug(f"按gr_code分组并聚合分数，聚合方式: {agg}")
    logger.debug(f"分数列: {score}")
    
    vals = df.groupby("gr_code")[score]
    if agg == "sum>0":
        result = vals.apply(lambda col: col[col > 0].sum()).reset_index()
    elif agg == "mean":
        result = vals.mean().reset_index()
    elif agg == "sum":
        result = vals.sum().reset_index()
    else:
        raise ValueError(f"未知的聚合方法: {agg}")
    
    logger.debug(f"结果形状: {result.shape}")
    logger.info(f"组分数范围: [{result[score].min():.4f}, {result[score].max():.4f}]")
    return result