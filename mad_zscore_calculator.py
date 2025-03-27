#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

def calculate_mad_zscore(values: pd.Series) -> pd.Series:
    """
    使用MAD计算Z分数
    
    参数:
        values: 输入数据序列
        
    返回:
        Z分数序列
    """
    median = values.median()
    mad = np.median(np.abs(values - median))
    return (values - median) / (1.4826 * mad)

def calculate_zrna(df: pd.DataFrame) -> pd.DataFrame:
    """
    计算ZRNA分数
    
    参数:
        df: 包含评估指标的DataFrame
        
    返回:
        添加了ZRNA分数的DataFrame
    """
    # 设置权重
    weights = {
        'TM': 1/3,
        'GDT': 1/3,
        'lDDT': 1/8,
        'INF_all': 1/8,
        'clash': 1/12
    }
    
    # 计算各个指标的Z分数
    for metric in weights.keys():
        df[f'{metric}_Zscore'] = calculate_mad_zscore(df[metric])
    
    # 计算综合Z分数
    zscore_components = []
    
    # 处理正向指标（越高越好）
    positive_metrics = ['TM', 'GDT', 'lDDT', 'INF_all']
    for metric in positive_metrics:
        component = weights[metric] * df[f'{metric}_Zscore']
        zscore_components.append(component)
    
    # 处理负向指标（越低越好）
    negative_metrics = ['clash']
    for metric in negative_metrics:
        component = -weights[metric] * df[f'{metric}_Zscore']
        zscore_components.append(component)
    
    # 计算最终的综合Z分数
    df['ZRNA'] = sum(zscore_components)
    
    return df

def main():
    # 设置输入输出路径
    input_file = "/Users/ooaahh/docs/ZscoreTest/all_RNA_AF3/R0250.csv"  # 替换为你的输入文件
    output_file = "/Users/ooaahh/docs/ZscoreTest/all_RNA_AF3/R0250_MAD_Zscore.csv"  # 替换为你的输出文件
    
    # 读取数据
    df = pd.read_csv(input_file)
    
    # 计算ZRNA
    result_df = calculate_zrna(df)
    
    # 保存结果
    result_df.to_csv(output_file, index=False)
    print(f"结果已保存到: {output_file}")

if __name__ == "__main__":
    main() 