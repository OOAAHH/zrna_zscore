data = readtable('R0250.csv', 'ReadVariableNames', true);

>> metrics = data{:, {'RMSD', 'INF_all', 'TM', 'GDT', 'clash', 'GDT', 'lDDT',}};

standardized_metrics = zscore(metrics);

new_column_names = strcat('Zscore_', metrics_columns);  % 将 'Zscore_' 前缀添加到每个列名

for i = 1:length(new_column_names)

    data.(new_column_names{i}) = standardized_metrics(:, i);  
end

weights = struct('TM', 1/3, 'GDT', 1/3, 'lDDT', 1/8, 'INF_all', 1/8, 'clash', 1/12);
TM_Z = data.Zscore_TM;
GDT_Z = data.Zscore_GDT;
lDDT_Z = data.Zscore_lDDT;
INF_all_Z = data.Zscore_INF_all;
clash_Z = data.Zscore_clash;
combined_Z = (weights.TM * TM_Z + weights.GDT * GDT_Z + ...
               weights.lDDT * lDDT_Z + weights.INF_all * INF_all_Z - ...
               weights.clash * clash_Z);
data.Combined_Zscore = combined_Z;