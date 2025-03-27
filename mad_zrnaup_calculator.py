#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path

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

def process_csv_file(file_path: str, output_dir: str):
    """
    处理单个csv文件
    
    参数:
        file_path: 输入文件路径
        output_dir: 输出目录
    """
    try:
        # 获取文件名（不含扩展名）
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        print(f"\n{'='*80}")
        print(f"开始处理文件: {file_name}")
        print(f"{'='*80}")
        
        # 读取csv文件
        print(f"正在读取文件: {file_path}")
        df = pd.read_csv(file_path)
        
        # 检查必要的列是否存在
        required_columns = ['TM', 'GDT', 'lDDT', 'INF_all', 'clash']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"错误: 文件 {file_name} 缺少必要的列: {missing_columns}")
            return
        
        # 计算ZRNA
        print(f"\n开始计算 {file_name} 的ZRNA分数")
        result_df = calculate_zrna(df)
        
        # 保存结果
        output_file = os.path.join(output_dir, f"{file_name}_MAD_ZRNA.csv")
        result_df.to_csv(output_file, index=False)
        print(f"\n结果已保存到: {output_file}")
        
        # 记录处理完成
        print(f"\n文件 {file_name} 处理完成")
        print(f"{'='*80}\n")
        
    except Exception as e:
        print(f"错误: 处理文件 {file_path} 时发生错误: {str(e)}")

def main():
    """
    主函数，处理指定目录下的所有csv文件
    """
    # 设置输入输出路径
    input_dir = "/Users/ooaahh/docs/ZscoreTest/logs/tmp"
    output_dir = "/Users/ooaahh/docs/ZscoreTest/logs/tmp_MAD_Zscore"
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    print(f"输出目录: {output_dir}")
    
    # 获取所有csv文件
    csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
    print(f"\n找到 {len(csv_files)} 个csv文件")
    for file in csv_files:
        print(f"  - {os.path.basename(file)}")
    
    # 处理每个文件
    for file_path in csv_files:
        process_csv_file(file_path, output_dir)
    
    print("\n" + "="*80)
    print("所有文件处理完成")
    print("="*80 + "\n")

if __name__ == "__main__":
    main() 