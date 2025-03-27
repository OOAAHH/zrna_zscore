function zscore_calculator()
    % Z分数计算程序
    % 用于计算RNA结构的综合Z分数
    
    % 设置输入输出路径
    input_dir = '/Users/ooaahh/docs/ZscoreTest/all_RNA_AF3';
    output_dir = '/Users/ooaahh/docs/ZscoreTest/zscore_results';
    log_dir = '/Users/ooaahh/docs/ZscoreTest/logs';
    
    % 创建输出和日志目录
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    if ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end
    
    % 设置日志文件
    log_file = fullfile(log_dir, 'zscore_calculation.log');
    diary(log_file);
    
    % 打印程序启动信息
    fprintf('='*80);
    fprintf('\nZ分数计算程序启动\n');
    fprintf('='*80);
    fprintf('\n');
    
    % 获取所有csv文件
    csv_files = dir(fullfile(input_dir, '*.csv'));
    fprintf('\n找到 %d 个csv文件\n', length(csv_files));
    for i = 1:length(csv_files)
        fprintf('  - %s\n', csv_files(i).name);
    end
    
    % 处理每个文件
    for i = 1:length(csv_files)
        process_csv_file(fullfile(input_dir, csv_files(i).name), output_dir);
    end
    
    % 打印程序结束信息
    fprintf('\n%s\n', '='*80);
    fprintf('所有文件处理完成\n');
    fprintf('%s\n', '='*80);
    
    % 关闭日志文件
    diary off;
end

function process_csv_file(file_path, output_dir)
    % 处理单个csv文件
    
    try
        % 获取文件名（不含扩展名）
        [~, file_name, ~] = fileparts(file_path);
        fprintf('\n%s\n', '='*80);
        fprintf('开始处理文件: %s\n', file_name);
        fprintf('%s\n', '='*80);
        
        % 读取csv文件
        fprintf('正在读取文件: %s\n', file_path);
        data = readtable(file_path);
        log_dataframe_info(data, '原始数据');
        
        % 检查必要的列是否存在
        required_columns = {'TM', 'GDT', 'lDDT', 'INF_all', 'clash'};
        missing_columns = setdiff(required_columns, data.Properties.VariableNames);
        if ~isempty(missing_columns)
            fprintf('错误: 文件 %s 缺少必要的列: %s\n', file_name, strjoin(missing_columns, ', '));
            return;
        end
        
        % 验证数据有效性
        for col = required_columns
            invalid_count = sum(isnan(data.(col{1})));
            if invalid_count > 0
                fprintf('警告: 列 %s 包含 %d 个无效值\n', col{1}, invalid_count);
            end
        end
        
        % 计算Z分数
        fprintf('\n开始计算 %s 的Z分数\n', file_name);
        result_data = calculate_comprehensive_zscore(data);
        
        % 保存结果
        output_file = fullfile(output_dir, sprintf('%s_zscore.csv', file_name));
        writetable(result_data, output_file, 'WriteRowNames', false);
        fprintf('\n结果已保存到: %s\n', output_file);
        log_dataframe_info(result_data, '计算结果');
        
        % 记录处理完成
        fprintf('\n文件 %s 处理完成\n', file_name);
        fprintf('%s\n', '='*80);
        
    catch e
        fprintf('错误: 处理文件 %s 时发生错误: %s\n', file_path, e.message);
    end
end

function log_dataframe_info(data, context)
    % 记录数据框的详细信息
    
    fprintf('\n%s %s %s\n', repmat('-', 1, 40), context, repmat('-', 1, 40));
    fprintf('数据形状: %d x %d\n', size(data));
    fprintf('列名列表:\n');
    for i = 1:length(data.Properties.VariableNames)
        fprintf('  - %s\n', data.Properties.VariableNames{i});
    end
    
    fprintf('\n数据预览:\n');
    disp(head(data, 5));
    
    fprintf('\n数据统计信息:\n');
    disp(summary(data));
    fprintf('%s\n', repmat('-', 1, 80));
end

function z_scores = calculate_zscore(values, metric_name)
    % 计算Z分数（标准分数）
    
    fprintf('\n计算 %s 的Z分数\n', metric_name);
    fprintf('输入数据统计信息:\n');
    fprintf('  数据长度: %d\n', length(values));
    fprintf('  数据范围: [%.4f, %.4f]\n', min(values), max(values));
    fprintf('  平均值: %.4f\n', mean(values));
    fprintf('  标准差: %.4f\n', std(values));
    
    % 检查数据是否包含NaN值
    nan_count = sum(isnan(values));
    if nan_count > 0
        fprintf('警告: 包含 %d 个NaN值\n', nan_count);
    end
    
    % 计算Z分数
    z_scores = (values - mean(values)) / std(values);
    
    fprintf('Z分数计算结果:\n');
    fprintf('  结果范围: [%.4f, %.4f]\n', min(z_scores), max(z_scores));
    fprintf('  平均值: %.4f\n', mean(z_scores));
    fprintf('  标准差: %.4f\n', std(z_scores));
end

function result_data = calculate_comprehensive_zscore(data)
    % 计算综合Z分数
    
    % 设置权重
    weights = struct(...
        'TM', 1/3, ...
        'GDT', 1/3, ...
        'lDDT', 1/8, ...
        'INF_all', 1/8, ...
        'clash', 1/12);
    
    fprintf('\n开始计算综合Z分数\n');
    fprintf('权重配置:\n');
    fields = fieldnames(weights);
    for i = 1:length(fields)
        fprintf('  %s: %.4f\n', fields{i}, weights.(fields{i}));
    end
    
    % 计算各个指标的Z分数
    for i = 1:length(fields)
        metric = fields{i};
        fprintf('\n处理指标: %s\n', metric);
        fprintf('  原始列名: %s\n', metric);
        fprintf('  目标列名: %s_Zscore\n', metric);
        data.(sprintf('%s_Zscore', metric)) = calculate_zscore(data.(metric), metric);
    end
    
    % 计算综合Z分数
    fprintf('\n计算综合Z分数\n');
    zscore_components = zeros(height(data), 1);
    
    % 处理正向指标（越高越好）
    positive_metrics = {'TM', 'GDT', 'lDDT', 'INF_all'};
    fprintf('\n处理正向指标（越高越好）:\n');
    for i = 1:length(positive_metrics)
        metric = positive_metrics{i};
        component = weights.(metric) * data.(sprintf('%s_Zscore', metric));
        zscore_components = zscore_components + component;
        fprintf('  %s:\n', metric);
        fprintf('    权重: %.4f\n', weights.(metric));
        fprintf('    贡献范围: [%.4f, %.4f]\n', min(component), max(component));
    end
    
    % 处理负向指标（越低越好）
    negative_metrics = {'clash'};
    fprintf('\n处理负向指标（越低越好）:\n');
    for i = 1:length(negative_metrics)
        metric = negative_metrics{i};
        component = -weights.(metric) * data.(sprintf('%s_Zscore', metric));
        zscore_components = zscore_components + component;
        fprintf('  %s:\n', metric);
        fprintf('    权重: %.4f\n', -weights.(metric));
        fprintf('    贡献范围: [%.4f, %.4f]\n', min(component), max(component));
    end
    
    % 添加综合Z分数到结果数据
    result_data = data;
    result_data.Zscore = zscore_components;
    
    fprintf('\n最终综合Z分数:\n');
    fprintf('  结果范围: [%.4f, %.4f]\n', min(zscore_components), max(zscore_components));
    fprintf('  平均值: %.4f\n', mean(zscore_components));
    fprintf('  标准差: %.4f\n', std(zscore_components));
end 