# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import logging
import os
import glob
from pathlib import Path
import shutil

def setup_logging(log_dir):
    """
    设置日志系统
    
    参数:
        log_dir (str): 日志文件存储目录
    """
    # 确保日志目录存在
    os.makedirs(log_dir, exist_ok=True)
    
    # 设置日志格式
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    # 设置日志处理器
    handlers = [
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(log_dir, 'zscore_calculation.log'))
    ]
    
    # 配置日志系统
    logging.basicConfig(
        level=logging.DEBUG,
        format=log_format,
        datefmt=date_format,
        handlers=handlers
    )
    
    logger = logging.getLogger(__name__)
    logger.info("="*80)
    logger.info("Z分数计算程序启动")
    logger.info("="*80)
    return logger

def log_dataframe_info(df: pd.DataFrame, logger, context: str):
    """
    记录DataFrame的详细信息
    
    参数:
        df: 要记录的DataFrame
        logger: 日志记录器
        context: 数据上下文说明
    """
    logger.info(f"\n{'-'*40} {context} {'-'*40}")
    logger.info(f"数据形状: {df.shape}")
    logger.info("列名列表:")
    for col in df.columns:
        logger.info(f"  - {col}")
    logger.info("\n数据预览:")
    logger.info(df.head())
    logger.info("\n数据统计信息:")
    logger.info(df.describe())
    logger.info(f"{'-'*80}\n")

def process_csv_file(file_path, output_dir, logger):
    """
    处理单个csv文件
    
    参数:
        file_path (str): csv文件路径
        output_dir (str): 输出目录
        logger: 日志记录器
    """
    try:
        # 获取文件名（不含扩展名）
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        logger.info(f"\n{'='*80}")
        logger.info(f"开始处理文件: {file_name}")
        logger.info(f"{'='*80}")
        
        # 读取csv文件
        logger.info(f"正在读取文件: {file_path}")
        df = pd.read_csv(file_path)
        log_dataframe_info(df, logger, "原始数据")
        
        # 检查必要的列是否存在
        required_columns = ['TM', 'GDT', 'lDDT', 'INF_all', 'clash']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.error(f"文件 {file_name} 缺少必要的列: {missing_columns}")
            return
        
        # 验证数据有效性
        for col in required_columns:
            invalid_count = df[col].isna().sum()
            if invalid_count > 0:
                logger.warning(f"列 {col} 包含 {invalid_count} 个无效值")
        
        # 计算Z分数
        logger.info(f"\n开始计算 {file_name} 的Z分数")
        result_df = calculate_comprehensive_zscore(df, logger)
        
        # 保存结果
        output_file = os.path.join(output_dir, f"{file_name}_zscore.csv")
        result_df.to_csv(output_file, index=False)
        logger.info(f"\n结果已保存到: {output_file}")
        log_dataframe_info(result_df, logger, "计算结果")
        
        # 记录处理完成
        logger.info(f"\n文件 {file_name} 处理完成")
        logger.info(f"{'='*80}\n")
        
    except Exception as e:
        logger.error(f"处理文件 {file_path} 时发生错误: {str(e)}", exc_info=True)

def calculate_zscore(values: pd.Series, logger, metric_name: str) -> pd.Series:
    """
    计算Z分数（标准分数）
    
    参数:
        values (pd.Series): 要标准化的Pandas系列数据
        logger: 日志记录器
        metric_name: 指标名称
        
    返回:
        pd.Series: 包含Z分数的Pandas系列
    """
    logger.info(f"\n计算 {metric_name} 的Z分数")
    logger.info(f"输入数据统计信息:")
    logger.info(f"  数据长度: {len(values)}")
    logger.info(f"  数据范围: [{values.min():.4f}, {values.max():.4f}]")
    logger.info(f"  平均值: {values.mean():.4f}")
    logger.info(f"  标准差: {values.std():.4f}")
    
    # 检查数据是否包含NaN值
    nan_count = values.isna().sum()
    if nan_count > 0:
        logger.warning(f"  包含 {nan_count} 个NaN值")
    
    # 计算均值和标准差
    mean_val = values.mean()
    std_val = values.std()
    
    # 标准化数值
    result = (values - mean_val) / std_val
    logger.info(f"Z分数计算结果:")
    logger.info(f"  结果范围: [{result.min():.4f}, {result.max():.4f}]")
    logger.info(f"  平均值: {result.mean():.4f}")
    logger.info(f"  标准差: {result.std():.4f}")
    
    return result

def calculate_comprehensive_zscore(
    df: pd.DataFrame,
    logger,
    zscore_column: str = 'Zscore',
    weights: dict = None
) -> pd.DataFrame:
    """
    计算综合Z分数
    
    参数:
        df (pd.DataFrame): 包含评估指标的DataFrame
        logger: 日志记录器
        zscore_column (str): 存储Z分数的列名
        weights (dict): 各指标的权重字典
        
    返回:
        pd.DataFrame: 添加了Z分数列的DataFrame
    """
    # 设置默认权重
    if weights is None:
        weights = {
            'TM': 1/3,
            'GDT': 1/3,
            'lDDT': 1/8,
            'INF_all': 1/8,
            'clash': 1/12
        }
    
    logger.info("\n开始计算综合Z分数")
    logger.info("权重配置:")
    for metric, weight in weights.items():
        logger.info(f"  {metric}: {weight:.4f}")
    
    # 计算各个指标的Z分数
    for metric in weights.keys():
        logger.info(f"\n处理指标: {metric}")
        logger.info(f"  原始列名: {metric}")
        logger.info(f"  目标列名: {metric}_Zscore")
        df[f'{metric}_Zscore'] = calculate_zscore(df[metric], logger, metric)
    
    # 计算综合Z分数
    logger.info("\n计算综合Z分数")
    zscore_components = []
    
    # 处理正向指标（越高越好）
    positive_metrics = ['TM', 'GDT', 'lDDT', 'INF_all']
    logger.info("\n处理正向指标（越高越好）:")
    for metric in positive_metrics:
        component = weights[metric] * df[f'{metric}_Zscore']
        zscore_components.append(component)
        logger.info(f"  {metric}:")
        logger.info(f"    权重: {weights[metric]:.4f}")
        logger.info(f"    贡献范围: [{component.min():.4f}, {component.max():.4f}]")
    
    # 处理负向指标（越低越好）
    negative_metrics = ['clash']
    logger.info("\n处理负向指标（越低越好）:")
    for metric in negative_metrics:
        component = -weights[metric] * df[f'{metric}_Zscore']
        zscore_components.append(component)
        logger.info(f"  {metric}:")
        logger.info(f"    权重: {-weights[metric]:.4f}")
        logger.info(f"    贡献范围: [{component.min():.4f}, {component.max():.4f}]")
    
    # 计算最终的综合Z分数
    df[zscore_column] = sum(zscore_components)
    logger.info(f"\n最终综合Z分数:")
    logger.info(f"  结果范围: [{df[zscore_column].min():.4f}, {df[zscore_column].max():.4f}]")
    logger.info(f"  平均值: {df[zscore_column].mean():.4f}")
    logger.info(f"  标准差: {df[zscore_column].std():.4f}")
    
    return df

def main():
    """
    主函数，处理指定目录下的所有csv文件
    """
    # 设置输入输出路径
    input_dir = "/Users/ooaahh/docs/ZscoreTest/all_RNA_AF3"
    output_dir = "/Users/ooaahh/docs/ZscoreTest/zscore_results"
    log_dir = "/Users/ooaahh/docs/ZscoreTest/logs"
    
    # 设置日志系统
    logger = setup_logging(log_dir)
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"输出目录: {output_dir}")
    
    # 获取所有csv文件
    csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
    logger.info(f"\n找到 {len(csv_files)} 个csv文件")
    for file in csv_files:
        logger.info(f"  - {os.path.basename(file)}")
    
    # 处理每个文件
    for file_path in csv_files:
        process_csv_file(file_path, output_dir, logger)
    
    logger.info("\n" + "="*80)
    logger.info("所有文件处理完成")
    logger.info("="*80 + "\n")

if __name__ == "__main__":
    main() 