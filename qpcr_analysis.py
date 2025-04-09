import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import time
from decimal import Decimal, ROUND_HALF_UP  # 添加decimal导入
import io  # 添加io模块用于处理文件流
import sys
import pkg_resources

def check_dependencies():
    """检查必要的依赖包是否已安装"""
    required_packages = {
        'streamlit': '1.24.0',
        'pandas': '1.5.0',
        'numpy': '1.21.0',
        'scipy': '1.9.0',
        'openpyxl': '3.0.10',
        'xlsxwriter': '3.0.9'
    }
    
    missing_packages = []
    for package, version in required_packages.items():
        try:
            pkg_resources.require(f"{package}>={version}")
        except (pkg_resources.DistributionNotFound, pkg_resources.VersionConflict):
            missing_packages.append(f"{package}>={version}")
    
    if missing_packages:
        st.error("缺少必要的依赖包，请运行以下命令安装：")
        st.code("pip install " + " ".join(missing_packages))
        sys.exit(1)

# 在程序启动时检查依赖
check_dependencies()

class QPCRAnalyzer:
    def calculate_standard_curve(self, concentrations, ct_values):
        """计算标准曲线"""
        slope, intercept, r_value, _, _ = stats.linregress(
            np.log10(concentrations), 
            ct_values
        )
        r_squared = r_value ** 2
        efficiency = (10 ** (-1/slope) - 1) * 100
        
        return slope, intercept, r_squared, efficiency
    
    def calculate_cv(self, values):
        """计算CV值"""
        return (np.std(values) / np.mean(values)) * 100

def clean_qpcr_data(df):
    """
    清洗QPCR数据，处理重复测量值
    - 所有基因共用一套Sample Name顺序
    """
    # 确保必要的列存在
    required_columns = ['Target Name', 'Sample Name', 'CT', 'Ct Mean', 'Quantity', 'Quantity Mean']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"缺少必要的列: {col}")
    
    # 选择需要的列
    df_cleaned = df[required_columns].copy()
    
    # 转换数值列的数据类型，保留Undetermined，将空值转为NA
    numeric_columns = ['CT', 'Ct Mean', 'Quantity', 'Quantity Mean']
    for col in numeric_columns:
        # 将空字符串转换为NA
        df_cleaned[col] = df_cleaned[col].replace('', 'NA')
        # 保持Undetermined不变
        mask_undetermined = df_cleaned[col] == 'Undetermined'
        # 将非Undetermined的值转换为Decimal并正确四舍五入到3位小数
        df_cleaned.loc[~mask_undetermined, col] = df_cleaned.loc[~mask_undetermined, col].apply(
            lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP)) 
            if pd.notnull(x) and x != 'NA' and x != 'Undetermined' 
            else x
        )
    
    # 获取所有唯一的Sample Name
    all_samples = df['Sample Name'].unique().tolist()
    
    # 通过文本输入设置Sample Name顺序
    st.write("设置所有基因共用的Sample Name顺序：")
    sample_order_input = st.text_area(
        "Sample Name顺序",
        height=200,
        help="每行输入一个Sample Name，未输入的样本将按原顺序排在后面"
    )
    
    # 处理输入的Sample Name顺序
    if sample_order_input:
        # 分割输入文本并去除空行和空白字符
        sample_order = [s.strip() for s in sample_order_input.split('\n') if s.strip()]
        
        # 验证输入的Sample Name是否有效
        invalid_samples = [s for s in sample_order if s not in all_samples]
        if invalid_samples:
            st.error(f"以下Sample Name不存在：{', '.join(invalid_samples)}")
            return None, None
        
        # 处理未被输入的样本
        remaining_samples = [s for s in all_samples if s not in sample_order]
        # 最终的样本顺序：用户输入的顺序 + 剩余样本的原始顺序
        final_sample_order = sample_order + remaining_samples
    else:
        # 如果没有输入，使用原始顺序
        final_sample_order = all_samples
    
    # 按Target Name分组处理数据
    target_names = df_cleaned['Target Name'].unique()
    target_dfs = {}
    
    for target in target_names:
        # 获取当前target的数据
        target_data = df_cleaned[df_cleaned['Target Name'] == target].copy()
        
        # 从原始数据中获取当前target的标准曲线参数
        target_params = df[df['Target Name'] == target]
        
        # 按Sample Name分组并重组数据
        grouped_data = []
        
        # 使用共同的Sample Name顺序处理数据
        for sample_name in final_sample_order:
            # 只处理当前target中存在的样本
            if sample_name in target_data['Sample Name'].values:
                group = target_data[target_data['Sample Name'] == sample_name]
                # 初始化行数据
                row = {'Sample Name': sample_name}
                
                # 处理CT值
                ct_values = group['CT'].tolist()
                for i, ct in enumerate(ct_values, 1):
                    if pd.isna(ct) or ct == '':
                        row[f'CT {i}'] = 'NA'
                    elif ct == 'Undetermined':
                        row[f'CT {i}'] = 'Undetermined'
                    else:
                        row[f'CT {i}'] = ct  # 保留原始值
                
                # 处理Ct Mean
                ct_mean = group['Ct Mean'].iloc[0] if not group['Ct Mean'].empty else 'NA'
                if pd.isna(ct_mean) or ct_mean == '':
                    row['Ct Mean'] = 'NA'
                elif ct_mean == 'Undetermined':
                    row['Ct Mean'] = 'Undetermined'
                else:
                    row['Ct Mean'] = ct_mean  # 保留原始值
                
                # 处理Quantity值
                quantity_values = group['Quantity'].tolist()
                for i, q in enumerate(quantity_values, 1):
                    if pd.isna(q) or q == '':
                        row[f'Quantity {i}'] = 'NA'
                    elif q == 'Undetermined':
                        row[f'Quantity {i}'] = 'Undetermined'
                    else:
                        row[f'Quantity {i}'] = q  # 保留原始值
                
                # 处理Quantity Mean
                quantity_mean = group['Quantity Mean'].iloc[0] if not group['Quantity Mean'].empty else 'NA'
                if pd.isna(quantity_mean) or quantity_mean == '':
                    row['Quantity Mean'] = 'NA'
                elif quantity_mean == 'Undetermined':
                    row['Quantity Mean'] = 'Undetermined'
                else:
                    row['Quantity Mean'] = quantity_mean  # 保留原始值
                
                # 从原始数据中获取当前target的R²值
                if 'R(superscript 2)' in df.columns:
                    r_squared = target_params[target_params['Sample Name'] == sample_name]['R(superscript 2)'].iloc[0]
                    if pd.isna(r_squared) or r_squared == '':
                        row['R²'] = 'NA'
                    else:
                        try:
                            # 直接格式化R²值为3位小数
                            row['R²'] = f"{float(r_squared):.3f}"
                        except (ValueError, TypeError):
                            row['R²'] = 'NA'
                else:
                    row['R²'] = 'NA'
                
                # 从原始数据中获取当前target的其他标准曲线参数
                # Y-Intercept
                if 'Y-Intercept' in df.columns:
                    y_intercept = target_params[target_params['Sample Name'] == sample_name]['Y-Intercept'].iloc[0]
                    if pd.isna(y_intercept) or y_intercept == '':
                        row['Y-Intercept'] = 'NA'
                    else:
                        try:
                            row['Y-Intercept'] = f"{float(y_intercept):.3f}"
                        except (ValueError, TypeError):
                            row['Y-Intercept'] = 'NA'
                else:
                    row['Y-Intercept'] = 'NA'
                
                # Slope
                if 'Slope' in df.columns:
                    slope = target_params[target_params['Sample Name'] == sample_name]['Slope'].iloc[0]
                    if pd.isna(slope) or slope == '':
                        row['Slope'] = 'NA'
                    else:
                        try:
                            row['Slope'] = f"{float(slope):.3f}"
                        except (ValueError, TypeError):
                            row['Slope'] = 'NA'
                else:
                    row['Slope'] = 'NA'
                
                # Efficiency
                if 'Efficiency' in df.columns:
                    efficiency = target_params[target_params['Sample Name'] == sample_name]['Efficiency'].iloc[0]
                    if pd.isna(efficiency) or efficiency == '':
                        row['Efficiency'] = 'NA'
                    else:
                        try:
                            row['Efficiency'] = f"{float(efficiency):.3f}"
                        except (ValueError, TypeError):
                            row['Efficiency'] = 'NA'
                else:
                    row['Efficiency'] = 'NA'
                
                grouped_data.append(row)
        
        # 创建DataFrame并设置列顺序
        result_df = pd.DataFrame(grouped_data)
        
        # 分离带数字的列和Mean列
        ct_numbered_cols = [col for col in result_df.columns if col.startswith('CT ')]
        quantity_numbered_cols = [col for col in result_df.columns if col.startswith('Quantity ') and 'Mean' not in col]
        
        # 按编号排序列名
        ct_numbered_cols.sort(key=lambda x: int(x.split()[1]))
        quantity_numbered_cols.sort(key=lambda x: int(x.split()[1]))
        
        # 设置最终列顺序，包括Quantity Mean
        col_order = (
            ['Sample Name'] + 
            ct_numbered_cols + 
            ['Ct Mean'] + 
            quantity_numbered_cols + 
            ['Quantity Mean'] +  # 添加Quantity Mean
            ['Y-Intercept', 'R²', 'Slope', 'Efficiency']
        )
        
        # 只选择存在的列
        col_order = [col for col in col_order if col in result_df.columns]
        result_df = result_df.reindex(columns=col_order)
        
        # 格式化显示（保留原始数据，仅在显示时格式化）
        for col in result_df.columns:
            if col.startswith(('CT', 'Quantity')):
                result_df[col] = result_df[col].apply(
                    lambda x: f"{float(x):.3f}" if isinstance(x, (int, float)) or (isinstance(x, str) and x.replace('.', '').isdigit()) else x
                )
        
        # 存储处理后的数据
        target_dfs[target] = result_df
    
    return target_dfs, df_cleaned

def calculate_cv_re(group_data):
    """计算CV和RE"""
    # 获取CT值列
    ct_cols = [col for col in group_data.columns if col.startswith('CT ') and col != 'CT Mean']
    
    if not ct_cols:
        return None, None
    
    # 转换为数值类型，排除'Undetermined'和'NA'
    ct_values = group_data[ct_cols].replace({'Undetermined': np.nan, 'NA': np.nan})
    ct_values = ct_values.astype(float)
    
    # 计算CV，保留1位小数
    cv = (ct_values.std(axis=1) / ct_values.mean(axis=1) * 100).round(1)
    
    # 计算RE (相对误差)，保留1位小数
    ct_mean = group_data['Ct Mean'].replace({'Undetermined': np.nan, 'NA': np.nan}).astype(float)
    re = ((ct_values.mean(axis=1) - ct_mean) / ct_mean * 100).round(1)
    
    return cv, re

def analyze_standard_curve(target_df, selected_target, selected_samples, cv_value_type, re_value_type, default_cv_threshold, default_re_threshold, custom_thresholds=None, theoretical_concentrations=None):
    """
    标准曲线分析
    """
    # 筛选选定的样本，并按输入顺序排序
    if selected_samples:
        # 创建一个字典，存储样本名称和其对应的顺序
        sample_order = {name: i for i, name in enumerate(selected_samples)}
        
        # 筛选数据并创建排序列
        target_df = target_df[target_df['Sample Name'].isin(selected_samples)].copy()
        target_df['sort_order'] = target_df['Sample Name'].map(sample_order)
        
        # 按照输入顺序排序
        target_df = target_df.sort_values('sort_order')
        
        # 删除排序列
        target_df = target_df.drop('sort_order', axis=1)
    
    # 创建结果DataFrame
    results = pd.DataFrame()
    results['Sample Name'] = target_df['Sample Name']
    
    # 1. 理论浓度
    if theoretical_concentrations:
        results['理论浓度'] = results['Sample Name'].map(theoretical_concentrations)
    else:
        results['理论浓度'] = np.nan
    
    # 2. 计算回算浓度
    def calculate_back_concentration_single(ct_value, slope, y_intercept):
        """计算单个CT值的回算浓度"""
        try:
            if ct_value in ['NA', 'Undetermined'] or pd.isna(ct_value):
                return 'NA'
            
            if slope is None or y_intercept is None:
                return 'NA'
            
            # 使用Decimal进行精确计算
            ct = Decimal(str(ct_value))
            slope_decimal = Decimal(str(slope))
            y_intercept_decimal = Decimal(str(y_intercept))
            
            # 使用标准曲线方程计算回算浓度：10^((Ct - b)/m)，其中b为y截距，m为斜率
            back_conc = Decimal('10') ** ((ct - y_intercept_decimal) / slope_decimal)
            
            # 保留3位小数
            return str(back_conc.quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
        except Exception as e:
            print(f"计算回算浓度时出错: {str(e)}")
            return 'NA'
    
    # 获取CT列
    ct_cols = [col for col in target_df.columns if col.startswith('CT ') and col != 'CT Mean']
    
    # 计算每个CT值对应的回算浓度
    for i, ct_col in enumerate(ct_cols, 1):
        results[f'回算浓度 {i}'] = target_df.apply(
            lambda row: calculate_back_concentration_single(
                row[ct_col],
                float(row['Slope']) if 'Slope' in row and row['Slope'] != 'NA' else None,
                float(row['Y-Intercept']) if 'Y-Intercept' in row and row['Y-Intercept'] != 'NA' else None
            ),
            axis=1
        )
    
    # 计算回算浓度Mean（3个回算浓度的平均值）
    def calculate_back_conc_mean(row):
        try:
            # 获取所有回算浓度值
            back_conc_values = []
            for i in range(1, 4):  # 假设有3个回算浓度
                value = row[f'回算浓度 {i}']
                if value not in ['NA', 'Undetermined'] and pd.notnull(value):
                    try:
                        back_conc_values.append(Decimal(str(value)))
                    except Exception as e:
                        print(f"回算浓度 {i} 转换错误: {value}, 错误: {str(e)}")
                        continue
            
            if not back_conc_values:
                print("没有有效的回算浓度值")
                return 'NA'
            
            # 使用Decimal进行精确计算平均值
            try:
                sum_value = sum(back_conc_values)
                count = Decimal(str(len(back_conc_values)))
                mean = sum_value / count
                
                # 保留3位小数
                result = str(mean.quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                print(f"回算浓度值: {back_conc_values}")
                print(f"计算得到的Mean: {result}")
                return result
            except Exception as e:
                print(f"计算平均值时出错: {str(e)}")
                return 'NA'
        except Exception as e:
            print(f"计算回算浓度Mean时出错: {str(e)}")
            return 'NA'
    
    results['回算浓度 Mean'] = results.apply(calculate_back_conc_mean, axis=1)
    
    # 3. 计算CV（使用CT值）
    cv = calculate_cv_by_type_standard(target_df)
    results['%CV'] = cv
    
    # 4. 计算RE（使用回算浓度Mean）
    re_series = pd.Series(index=results.index)
    for idx, row in results.iterrows():
        try:
            sample_name = row['Sample Name']
            theo_conc = theoretical_concentrations.get(sample_name, np.nan)
            
            if pd.isna(theo_conc):
                re_series[idx] = np.nan
                continue
            
            # 使用回算浓度Mean
            back_conc_mean = row['回算浓度 Mean']
            
            if back_conc_mean in ['NA', 'Undetermined'] or pd.isna(back_conc_mean):
                re_series[idx] = np.nan
                continue
            
            # 确保数值转换正确
            try:
                theo_conc_decimal = Decimal(str(theo_conc))
                back_conc_decimal = Decimal(str(back_conc_mean))
            except Exception as e:
                print(f"数值转换错误 - 样本 {sample_name}: {str(e)}")
                print(f"理论浓度: {theo_conc}, 回算浓度Mean: {back_conc_mean}")
                re_series[idx] = np.nan
                continue
            
            if theo_conc_decimal == 0:
                re_series[idx] = np.nan
                continue
            
            # 计算RE: (回算浓度Mean - 理论浓度) / 理论浓度 * 100
            try:
                diff = back_conc_decimal - theo_conc_decimal
                ratio = diff / theo_conc_decimal
                re = (ratio * Decimal('100')).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)
                re_series[idx] = str(re)
            except Exception as e:
                print(f"RE计算错误 - 样本 {sample_name}: {str(e)}")
                print(f"理论浓度: {theo_conc_decimal}, 回算浓度Mean: {back_conc_decimal}")
                re_series[idx] = np.nan
            
        except Exception as e:
            print(f"计算RE时出错 - 样本 {sample_name}: {str(e)}")
            re_series[idx] = np.nan
    
    results['%RE'] = re_series
    
    # 5. 判定结果
    def check_acceptance(row):
        try:
            if pd.isna(row['%CV']) or pd.isna(row['%RE']):
                return 'NA'
            
            # 使用Decimal进行精确比较
            cv_value = Decimal(str(row['%CV'])) if isinstance(row['%CV'], str) else Decimal(str(row['%CV']))
            re_value = Decimal(str(row['%RE'])) if isinstance(row['%RE'], str) else Decimal(str(row['%RE']))
            
            if cv_value <= Decimal(str(default_cv_threshold)) and abs(re_value) <= Decimal(str(default_re_threshold)):
                return '符合标准'
            return '不符合标准'
        except (ValueError, TypeError) as e:
            print(f"判定结果时出错: {str(e)}")
            return 'NA'
    
    results['判定结果'] = results.apply(check_acceptance, axis=1)
    
    # 格式化数值显示（保持Decimal字符串格式）
    for col in ['%CV', '%RE']:
        results[col] = results[col].apply(
            lambda x: str(Decimal(str(x)).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)) 
            if pd.notnull(x) and x != 'NA' and x != 'Undetermined' 
            else x
        )
    
    # 格式化理论浓度和回算浓度显示（保持Decimal字符串格式）
    for col in results.columns:
        if col == '理论浓度' or col.startswith('回算浓度'):
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else 'NA'
            )
    
    return results, None

def calculate_cv_by_type(target_df, value_type):
    """根据选择的计算方式计算CV"""
    if value_type == "CT值":
        # 获取CT列
        ct_cols = [col for col in target_df.columns if col.startswith('CT ') and col != 'CT Mean']
        
        def calculate_cv_row(row):
            try:
                # 获取有效的CT值和CT Mean
                ct_values = []
                for col in ct_cols:
                    value = row[col]
                    if value not in ['NA', 'Undetermined'] and pd.notnull(value):
                        ct_values.append(Decimal(str(value)))
                
                ct_mean = row['Ct Mean']
                if ct_mean in ['NA', 'Undetermined'] or pd.isna(ct_mean):
                    return np.nan
                
                ct_mean = Decimal(str(ct_mean))
                if ct_mean == 0:
                    return np.nan
                
                if not ct_values:
                    return np.nan
                
                # 计算标准差
                squared_diff_sum = sum((x - ct_mean) ** 2 for x in ct_values)
                n = len(ct_values)
                if n <= 1:
                    return np.nan
                    
                std = (squared_diff_sum / Decimal(str(n - 1))).sqrt()
                
                # 计算CV并保留1位小数
                cv = (std / ct_mean * Decimal('100')).quantize(
                    Decimal('0.1'), 
                    rounding=ROUND_HALF_UP
                )
                return str(cv)  # 返回字符串形式的Decimal值
            except Exception as e:
                print(f"计算CV时出错: {str(e)}")
                return np.nan
        
        cv = target_df.apply(calculate_cv_row, axis=1)
        return cv
    
    return pd.Series(index=target_df.index)

def calculate_re_by_type_standard(target_df, theoretical_concentrations=None):
    """标准曲线分析的RE计算（使用回算浓度Mean）"""
    if theoretical_concentrations is None:
        return pd.Series(index=target_df.index)
    
    re_series = pd.Series(index=target_df.index)
    
    # 计算RE
    for idx, row in target_df.iterrows():
        try:
            sample_name = row['Sample Name']
            theo_conc = theoretical_concentrations.get(sample_name, np.nan)
            
            if pd.isna(theo_conc):
                re_series[idx] = np.nan
                continue
            
            # 使用回算浓度Mean
            back_conc_mean = row['回算浓度 Mean']
            
            if back_conc_mean in ['NA', 'Undetermined'] or pd.isna(back_conc_mean):
                re_series[idx] = np.nan
                continue
            
            theo_conc_decimal = Decimal(str(theo_conc))
            back_conc_decimal = Decimal(str(back_conc_mean))
            
            if theo_conc_decimal == 0:
                re_series[idx] = np.nan
                continue
            
            # 计算RE并保留1位小数
            re = ((back_conc_decimal - theo_conc_decimal) / theo_conc_decimal * Decimal('100')).quantize(
                Decimal('0.1'),
                rounding=ROUND_HALF_UP
            )
            re_series[idx] = str(re)  # 返回字符串形式的Decimal值
            
        except Exception as e:
            print(f"计算RE时出错: {str(e)}")
            re_series[idx] = np.nan
    
    return re_series

def calculate_cv_by_type_standard(target_df):
    """标准曲线分析的CV计算（使用CT值）"""
    # 获取CT列
    ct_cols = [col for col in target_df.columns if col.startswith('CT ') and col != 'CT Mean']
    
    def calculate_cv_row(row):
        try:
            # 获取有效的CT值和CT Mean
            ct_values = []
            for col in ct_cols:
                value = row[col]
                if value not in ['NA', 'Undetermined'] and pd.notnull(value):
                    ct_values.append(Decimal(str(value)))
            
            ct_mean = row['Ct Mean']
            if ct_mean in ['NA', 'Undetermined'] or pd.isna(ct_mean):
                return np.nan
            
            ct_mean = Decimal(str(ct_mean))
            if ct_mean == 0:
                return np.nan
            
            if not ct_values:
                return np.nan
            
            # 计算标准差
            squared_diff_sum = sum((x - ct_mean) ** 2 for x in ct_values)
            n = len(ct_values)
            if n <= 1:
                return np.nan
                
            std = (squared_diff_sum / Decimal(str(n - 1))).sqrt()
            
            # 计算CV并保留1位小数
            cv = (std / ct_mean * Decimal('100')).quantize(
                Decimal('0.1'), 
                rounding=ROUND_HALF_UP
            )
            return str(cv)  # 返回字符串形式的Decimal值
        except Exception as e:
            print(f"计算CV时出错: {str(e)}")
            return np.nan
    
    return target_df.apply(calculate_cv_row, axis=1)

def analyze_control_samples(target_df, selected_samples, value_type, threshold_type, threshold=None, min_threshold=None, max_threshold=None):
    """
    对照样品分析
    - 显示时：空值显示为NA
    - 计算时：CT/CT Mean空值为40，Quantity/Quantity Mean空值为0
    """
    # 筛选选定的样本
    if selected_samples:
        target_df = target_df[target_df['Sample Name'].isin(selected_samples)].copy()
    
    # 创建结果DataFrame的基础数据
    results_data = {
        'Sample Name': target_df['Sample Name']
    }
    
    # 根据值类型选择要分析的列
    if value_type == "CT值":
        # 获取所有CT列
        value_cols = [col for col in target_df.columns if col.startswith('CT ')]
        mean_col = 'Ct Mean'
        default_value = 40  # CT值的默认值
    else:  # Quantity
        # 获取所有Quantity列
        value_cols = [col for col in target_df.columns if col.startswith('Quantity ')]
        mean_col = 'Quantity Mean'
        default_value = 0  # Quantity的默认值
    
    # 添加标准范围信息
    if threshold_type == "范围":
        results_data['标准范围'] = f"{min_threshold:.3f} - {max_threshold:.3f}"
    elif threshold_type == "小于":
        results_data['标准范围'] = f"< {threshold:.3f}"
    elif threshold_type == "小于等于":
        results_data['标准范围'] = f"≤ {threshold:.3f}"
    elif threshold_type == "大于":
        results_data['标准范围'] = f"> {threshold:.3f}"
    else:  # 大于等于
        results_data['标准范围'] = f"≥ {threshold:.3f}"
    
    # 添加所有值列（保持显示格式）
    for col in value_cols:
        results_data[col] = target_df[col]
    
    results = pd.DataFrame(results_data)
    
    # 添加判断结果
    def check_acceptance(row):
        try:
            # 检查所有值
            all_values = []
            for col in value_cols:
                value = row[col]
                # 计算时将NA和Undetermined转换为默认值
                if value in ['NA', 'Undetermined'] or pd.isna(value):
                    calc_value = default_value
                else:
                    calc_value = float(value)
                all_values.append(calc_value)
            
            # 根据标准类型进行判断
            if threshold_type == "范围":
                # 检查是否所有值都在范围内
                if all(min_threshold <= q <= max_threshold for q in all_values):
                    return "符合标准"
                else:
                    # 创建每个值的判定结果
                    value_results = []
                    for i, q in enumerate(all_values, 1):
                        if min_threshold <= q <= max_threshold:
                            value_results.append(f"{value_type}{i}: 符合")
                        else:
                            value_results.append(f"{value_type}{i}: 不符合")
                    return f"不符合标准 ({'; '.join(value_results)})"
            
            elif threshold_type == "小于":
                # 检查是否所有值都小于阈值
                if all(q < threshold for q in all_values):
                    return "符合标准"
                else:
                    # 创建每个值的判定结果
                    value_results = []
                    for i, q in enumerate(all_values, 1):
                        if q < threshold:
                            value_results.append(f"{value_type}{i}: 符合")
                        else:
                            value_results.append(f"{value_type}{i}: 不符合")
                    return f"不符合标准 ({'; '.join(value_results)})"
            
            elif threshold_type == "小于等于":
                # 检查是否所有值都小于等于阈值
                if all(q <= threshold for q in all_values):
                    return "符合标准"
                else:
                    # 创建每个值的判定结果
                    value_results = []
                    for i, q in enumerate(all_values, 1):
                        if q <= threshold:
                            value_results.append(f"{value_type}{i}: 符合")
                        else:
                            value_results.append(f"{value_type}{i}: 不符合")
                    return f"不符合标准 ({'; '.join(value_results)})"
            
            elif threshold_type == "大于":
                # 检查是否所有值都大于阈值
                if all(q > threshold for q in all_values):
                    return "符合标准"
                else:
                    # 创建每个值的判定结果
                    value_results = []
                    for i, q in enumerate(all_values, 1):
                        if q > threshold:
                            value_results.append(f"{value_type}{i}: 符合")
                        else:
                            value_results.append(f"{value_type}{i}: 不符合")
                    return f"不符合标准 ({'; '.join(value_results)})"
            
            else:  # 大于等于
                # 检查是否所有值都大于等于阈值
                if all(q >= threshold for q in all_values):
                    return "符合标准"
                else:
                    # 创建每个值的判定结果
                    value_results = []
                    for i, q in enumerate(all_values, 1):
                        if q >= threshold:
                            value_results.append(f"{value_type}{i}: 符合")
                        else:
                            value_results.append(f"{value_type}{i}: 不符合")
                    return f"不符合标准 ({'; '.join(value_results)})"
                
        except (ValueError, TypeError):
            return 'NA'
    
    results['判定结果'] = results.apply(check_acceptance, axis=1)
    
    return results

def calculate_re_by_type_qc(target_df, theoretical_concentrations=None):
    """质控样品分析的RE计算（使用Quantity Mean）"""
    if theoretical_concentrations is None:
        return pd.Series(index=target_df.index)
    
    re_series = pd.Series(index=target_df.index)
    
    for idx, row in target_df.iterrows():
        try:
            sample_name = row['Sample Name']
            theo_conc = theoretical_concentrations.get(sample_name, np.nan)
            
            if pd.isna(theo_conc):
                re_series[idx] = np.nan
                continue
            
            # 使用Quantity Mean
            quantity_mean = row['Quantity Mean']
            
            if quantity_mean in ['NA', 'Undetermined'] or pd.isna(quantity_mean):
                re_series[idx] = np.nan
                continue
            
            theo_conc_decimal = Decimal(str(theo_conc))
            quantity_mean_decimal = Decimal(str(quantity_mean))
            
            if theo_conc_decimal == 0:
                re_series[idx] = np.nan
                continue
            
            # 计算RE并保留1位小数
            re = ((quantity_mean_decimal - theo_conc_decimal) / theo_conc_decimal * Decimal('100')).quantize(
                Decimal('0.1'),
                rounding=ROUND_HALF_UP
            )
            re_series[idx] = str(re)  # 返回字符串形式的Decimal值
            
        except Exception as e:
            print(f"计算RE时出错: {str(e)}")
            re_series[idx] = np.nan
    
    return re_series

def calculate_cv_by_type_qc(target_df):
    """质控样品分析的CV计算（使用Quantity值）"""
    # 获取Quantity列
    quantity_cols = [col for col in target_df.columns if col.startswith('Quantity ')]
    
    def calculate_cv_row(row):
        try:
            # 获取有效的Quantity值
            quantity_values = []
            for col in quantity_cols:
                value = row[col]
                if value not in ['NA', 'Undetermined'] and pd.notnull(value):
                    quantity_values.append(Decimal(str(value)))
            
            if not quantity_values:
                return np.nan
            
            # 计算标准差
            squared_diff_sum = sum((x - sum(quantity_values)/Decimal(str(len(quantity_values)))) ** 2 for x in quantity_values)
            n = len(quantity_values)
            if n <= 1:
                return np.nan
                
            std = (squared_diff_sum / Decimal(str(n - 1))).sqrt()
            mean = sum(quantity_values) / Decimal(str(n))
            
            if mean == 0:
                return np.nan
            
            # 计算CV并保留1位小数
            cv = (std / mean * Decimal('100')).quantize(
                Decimal('0.1'), 
                rounding=ROUND_HALF_UP
            )
            return str(cv)  # 返回字符串形式的Decimal值
        except Exception as e:
            print(f"计算CV时出错: {str(e)}")
            return np.nan
    
    return target_df.apply(calculate_cv_row, axis=1)

def analyze_qc_samples(target_df, selected_target, selected_samples, cv_value_type, re_value_type, default_cv_threshold, default_re_threshold, custom_thresholds=None, theoretical_concentrations=None):
    """质控样品分析"""
    # 筛选选定的样本，并按输入顺序排序
    if selected_samples:
        # 创建一个字典，存储样本名称和其对应的顺序
        sample_order = {name: i for i, name in enumerate(selected_samples)}
        
        # 筛选数据并创建排序列
        target_df = target_df[target_df['Sample Name'].isin(selected_samples)].copy()
        target_df['sort_order'] = target_df['Sample Name'].map(sample_order)
        
        # 按照输入顺序排序
        target_df = target_df.sort_values('sort_order')
        
        # 删除排序列
        target_df = target_df.drop('sort_order', axis=1)
    
    # 创建结果DataFrame
    results = pd.DataFrame()
    results['Sample Name'] = target_df['Sample Name']
    
    # 添加理论浓度（放在Sample Name后面）
    if theoretical_concentrations:
        results['理论浓度'] = results['Sample Name'].map(theoretical_concentrations)
    else:
        results['理论浓度'] = np.nan
    
    # 添加Quantity列和Quantity Mean列
    quantity_cols = [col for col in target_df.columns if col.startswith('Quantity ')]
    for col in quantity_cols:
        results[col] = target_df[col]
    results['Quantity Mean'] = target_df['Quantity Mean']
    
    # 计算CV（使用Quantity）
    cv = calculate_cv_by_type_qc(target_df)
    results['%CV'] = cv
    
    # 计算RE（使用Quantity Mean）
    re = calculate_re_by_type_qc(target_df, theoretical_concentrations)
    results['%RE'] = re
    
    # 判定结果
    def check_acceptance(row):
        try:
            if pd.isna(row['%CV']) or pd.isna(row['%RE']):
                return 'NA'
            
            # 使用Decimal进行精确比较
            cv_value = Decimal(str(row['%CV'])) if isinstance(row['%CV'], str) else Decimal(str(row['%CV']))
            re_value = Decimal(str(row['%RE'])) if isinstance(row['%RE'], str) else Decimal(str(row['%RE']))
            
            if cv_value <= Decimal(str(default_cv_threshold)) and abs(re_value) <= Decimal(str(default_re_threshold)):
                return '符合标准'
            return '不符合标准'
        except (ValueError, TypeError) as e:
            print(f"判定结果时出错: {str(e)}")
            return 'NA'
    
    results['判定结果'] = results.apply(check_acceptance, axis=1)
    
    # 格式化数值显示（保持Decimal字符串格式）
    for col in ['%CV', '%RE']:
        results[col] = results[col].apply(
            lambda x: str(Decimal(str(x)).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)) 
            if pd.notnull(x) and x != 'NA' and x != 'Undetermined' 
            else x
        )
    
    # 格式化理论浓度、Quantity和Quantity Mean显示（保持Decimal字符串格式）
    for col in ['理论浓度'] + quantity_cols + ['Quantity Mean']:
        results[col] = results[col].apply(
            lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
            if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
            else 'NA'
        )
    
    return results, None

def analyze_test_samples(target_df, selected_target, selected_samples, value_type, threshold_type=None, threshold=None, min_threshold=None, max_threshold=None, cv_threshold=None):
    """待测样品分析"""
    # 筛选选定的样本，并按输入顺序排序
    if selected_samples:
        # 创建一个字典，存储样本名称和其对应的顺序
        sample_order = {name: i for i, name in enumerate(selected_samples)}
        
        # 筛选数据并创建排序列
        target_df = target_df[target_df['Sample Name'].isin(selected_samples)].copy()
        target_df['sort_order'] = target_df['Sample Name'].map(sample_order)
        
        # 按照输入顺序排序
        target_df = target_df.sort_values('sort_order')
        
        # 删除排序列
        target_df = target_df.drop('sort_order', axis=1)
    
    # 创建结果DataFrame
    results = pd.DataFrame()
    results['Sample Name'] = target_df['Sample Name']
    
    # 根据选择的值类型添加相应的列
    if value_type in ["CT值", "CT值和Quantity"]:
        # 添加CT列和CT Mean列
        ct_cols = [col for col in target_df.columns if col.startswith('CT ')]
        for col in ct_cols:
            results[col] = target_df[col]
        results['Ct Mean'] = target_df['Ct Mean']
        
        # 计算CT的CV
        cv_ct = calculate_cv_by_type_standard(target_df)
        results['CT %CV'] = cv_ct
    
    if value_type in ["Quantity", "CT值和Quantity"]:
        # 添加Quantity列和Quantity Mean列
        quantity_cols = [col for col in target_df.columns if col.startswith('Quantity ')]
        for col in quantity_cols:
            results[col] = target_df[col]
        results['Quantity Mean'] = target_df['Quantity Mean']
        
        # 计算Quantity的CV
        cv_quantity = calculate_cv_by_type_qc(target_df)
        results['Quantity %CV'] = cv_quantity
    
    # 添加标准范围信息
    if threshold_type:
        if threshold_type == "范围":
            results['标准范围'] = f"{min_threshold:.3f} - {max_threshold:.3f}"
        elif threshold_type == "小于":
            results['标准范围'] = f"< {threshold:.3f}"
        elif threshold_type == "小于等于":
            results['标准范围'] = f"≤ {threshold:.3f}"
        elif threshold_type == "大于":
            results['标准范围'] = f"> {threshold:.3f}"
        else:  # 大于等于
            results['标准范围'] = f"≥ {threshold:.3f}"
    
    # 添加CV接受标准信息
    if cv_threshold is not None:
        results['CV接受标准'] = f"≤ {cv_threshold:.1f}%"
        
        # 添加CV判定结果
        def check_cv_acceptance(row):
            try:
                if value_type == "CT值":
                    cv_value = row['CT %CV']
                elif value_type == "Quantity":
                    cv_value = row['Quantity %CV']
                else:  # CT值和Quantity
                    cv_ct = row['CT %CV']
                    cv_quantity = row['Quantity %CV']
                    
                    if cv_ct == 'NA' and cv_quantity == 'NA':
                        return 'NA'
                    elif cv_ct == 'NA':
                        cv_value = cv_quantity
                    elif cv_quantity == 'NA':
                        cv_value = cv_ct
                    else:
                        cv_value = max(float(cv_ct), float(cv_quantity))
                
                if cv_value == 'NA':
                    return 'NA'
                
                cv_decimal = Decimal(str(cv_value))
                threshold_decimal = Decimal(str(cv_threshold))
                
                if cv_decimal <= threshold_decimal:
                    return "符合标准"
                return "不符合标准"
                
            except (ValueError, TypeError):
                return 'NA'
        
        results['CV判定结果'] = results.apply(check_cv_acceptance, axis=1)
    
    # 添加判定结果
    if threshold_type:
        def check_acceptance(row):
            try:
                if value_type in ["CT值", "CT值和Quantity"]:
                    # 检查CT值
                    ct_values = []
                    for col in ct_cols:
                        value = row[col]
                        if value not in ['NA', 'Undetermined'] and pd.notnull(value):
                            ct_values.append(float(value))
                    
                    if ct_values:
                        if threshold_type == "范围":
                            ct_valid = all(min_threshold <= q <= max_threshold for q in ct_values)
                        elif threshold_type == "小于":
                            ct_valid = all(q < threshold for q in ct_values)
                        elif threshold_type == "小于等于":
                            ct_valid = all(q <= threshold for q in ct_values)
                        elif threshold_type == "大于":
                            ct_valid = all(q > threshold for q in ct_values)
                        else:  # 大于等于
                            ct_valid = all(q >= threshold for q in ct_values)
                    else:
                        ct_valid = False
                
                if value_type in ["Quantity", "CT值和Quantity"]:
                    # 检查Quantity值
                    quantity_values = []
                    for col in quantity_cols:
                        value = row[col]
                        if value not in ['NA', 'Undetermined'] and pd.notnull(value):
                            quantity_values.append(float(value))
                    
                    if quantity_values:
                        if threshold_type == "范围":
                            quantity_valid = all(min_threshold <= q <= max_threshold for q in quantity_values)
                        elif threshold_type == "小于":
                            quantity_valid = all(q < threshold for q in quantity_values)
                        elif threshold_type == "小于等于":
                            quantity_valid = all(q <= threshold for q in quantity_values)
                        elif threshold_type == "大于":
                            quantity_valid = all(q > threshold for q in quantity_values)
                        else:  # 大于等于
                            quantity_valid = all(q >= threshold for q in quantity_values)
                    else:
                        quantity_valid = False
                
                # 根据选择的值类型返回判定结果
                if value_type == "CT值":
                    return "符合标准" if ct_valid else "不符合标准"
                elif value_type == "Quantity":
                    return "符合标准" if quantity_valid else "不符合标准"
                else:  # CT值和Quantity
                    if ct_valid and quantity_valid:
                        return "符合标准"
                    elif not ct_valid and not quantity_valid:
                        return "不符合标准"
                    else:
                        invalid_types = []
                        if not ct_valid:
                            invalid_types.append("CT值")
                        if not quantity_valid:
                            invalid_types.append("Quantity")
                        return f"部分不符合标准 ({', '.join(invalid_types)})"
                
            except (ValueError, TypeError):
                return 'NA'
        
        results['判定结果'] = results.apply(check_acceptance, axis=1)
    
    # 格式化数值显示（保持Decimal字符串格式）
    for col in results.columns:
        if col.startswith(('CT', 'Quantity')):
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else 'NA'
            )
        elif col.endswith('%CV'):
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else x
            )
    
    # 确保所有列的数据格式一致
    for col in results.columns:
        if col not in ['Sample Name', '标准范围', '判定结果', '阴阳性判定']:
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else 'NA'
            )
    
    return results, None

def generate_report(target_dfs, selected_target, selected_samples=None, report_type="VCN", determination_type=None, determination_rule=None, threshold_value=None, min_value=None, max_value=None):
    """
    生成汇报数据
    report_type: 汇报类型，可选 "VCN", "RCL", "DNA中拷贝数"
    determination_type: 判定方式，可选 "CT值", "Quantity值"
    determination_rule: 判定规则，可选 "小于", "小于等于", "大于", "大于等于", "范围"
    threshold_value: 阈值（当判定规则不是"范围"时使用）
    min_value: 最小值（当判定规则为"范围"时使用）
    max_value: 最大值（当判定规则为"范围"时使用）
    """
    if selected_target not in target_dfs:
        return None
    
    target_df = target_dfs[selected_target]
    
    # 筛选选定的样本
    if selected_samples:
        target_df = target_df[target_df['Sample Name'].isin(selected_samples)].copy()
    
    # 创建结果DataFrame
    results = pd.DataFrame()
    results['Sample Name'] = target_df['Sample Name']
    
    # 根据汇报类型生成不同的报告
    if report_type == "VCN":
        # VCN汇报：显示CT值和Quantity值
        ct_cols = [col for col in target_df.columns if col.startswith('CT ')]
        quantity_cols = [col for col in target_df.columns if col.startswith('Quantity ')]
        
        # 添加CT列
        for col in ct_cols:
            results[col] = target_df[col]
        results['Ct Mean'] = target_df['Ct Mean']
        
        # 添加Quantity列
        for col in quantity_cols:
            results[col] = target_df[col]
        results['Quantity Mean'] = target_df['Quantity Mean']
        
        # 计算CV
        cv = calculate_cv_by_type_standard(target_df)
        results['%CV'] = cv
    
    elif report_type == "RCL":
        # RCL汇报：只显示Sample Name、CT、Ct Mean和阴阳性判定
        ct_cols = [col for col in target_df.columns if col.startswith('CT ')]
        
        # 添加CT列
        for col in ct_cols:
            results[col] = target_df[col]
        results['Ct Mean'] = target_df['Ct Mean']
        
        # 添加阴阳性判定（使用Ct Mean判断）
        def check_determination(row):
            try:
                ct_mean = row['Ct Mean']
                
                # 处理无效值
                if ct_mean in ['NA', 'Undetermined'] or pd.isna(ct_mean):
                    return '无效'
                
                ct_mean = float(ct_mean)
                
                # 根据判定规则进行判断
                if determination_rule == "范围":
                    if min_value <= ct_mean <= max_value:
                        return '阳性'
                    return '阴性'
                elif determination_rule == "小于":
                    if ct_mean < threshold_value:
                        return '阳性'
                    return '阴性'
                elif determination_rule == "小于等于":
                    if ct_mean <= threshold_value:
                        return '阳性'
                    return '阴性'
                elif determination_rule == "大于":
                    if ct_mean > threshold_value:
                        return '阳性'
                    return '阴性'
                else:  # 大于等于
                    if ct_mean >= threshold_value:
                        return '阳性'
                    return '阴性'
                    
            except (ValueError, TypeError):
                return '无效'
        
        results['阴阳性判定'] = results.apply(check_determination, axis=1)
    
    else:  # DNA中拷贝数
        # DNA中拷贝数汇报：显示Quantity值和拷贝数
        quantity_cols = [col for col in target_df.columns if col.startswith('Quantity ')]
        
        # 添加Quantity列
        for col in quantity_cols:
            results[col] = target_df[col]
        results['Quantity Mean'] = target_df['Quantity Mean']
        
        # 计算拷贝数（假设每个细胞有2个拷贝）
        results['拷贝数'] = results.apply(
            lambda row: float(row['Quantity Mean']) * 2 
            if row['Quantity Mean'] not in ['NA', 'Undetermined'] 
            else 'NA',
            axis=1
        )
        
        # 计算CV
        cv = calculate_cv_by_type_qc(target_df)
        results['%CV'] = cv
    
    # 格式化数值显示
    for col in results.columns:
        if col.startswith(('CT', 'Quantity')):
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else 'NA'
            )
        elif col.endswith('%CV'):
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else x
            )
        elif col in ['相对表达量', '拷贝数']:
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else 'NA'
            )
    
    # 确保所有列的数据格式一致
    for col in results.columns:
        if col not in ['Sample Name', '标准范围', '判定结果', '阴阳性判定']:
            results[col] = results[col].apply(
                lambda x: str(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                if pd.notnull(x) and x != 'NA' and x != 'Undetermined'
                else 'NA'
            )
    
    return results

def process_module_data(result, analysis_name, is_last_module=False):
    """处理单个模块的数据，返回处理后的DataFrame"""
    all_results = pd.DataFrame()
    
    # 添加模块标题（只添加标题行）
    module_title = pd.DataFrame([
        [analysis_name]  # 直接使用模块名称作为标题
    ], columns=['Sample Name'])
    all_results = pd.concat([all_results, module_title], ignore_index=True)
    
    if analysis_name == '对照样品分析':
        # 处理对照样品分析的多组结果
        for i, group_result in enumerate(result):
            # 删除全是NA的列
            valid_cols = ['Sample Name']  # 始终保留Sample Name列
            for col in group_result.columns:
                if col != 'Sample Name':
                    # 检查列中是否所有值都是NA或Undetermined
                    values = group_result[col].astype(str)
                    if not all(values.isin(['NA', 'Undetermined', 'nan', 'None', ''])):
                        valid_cols.append(col)
            group_result = group_result[valid_cols]
            
            # 添加组标题（只添加标题行）
            group_title = pd.DataFrame([
                [f"第{i+1}组对照分析结果"]  # 组标题
            ], columns=['Sample Name'])
            all_results = pd.concat([all_results, group_title], ignore_index=True)
            
            # 添加数据（包括表头）
            group_data = pd.DataFrame(columns=group_result.columns)
            group_data.loc[0] = group_result.columns  # 添加表头
            
            # 处理数值数据
            processed_group_result = group_result.copy()
            for col in processed_group_result.columns:
                if col != 'Sample Name':
                    processed_group_result[col] = processed_group_result[col].apply(
                        lambda x: str(x) if pd.isna(x) or not isinstance(x, (int, float)) or x in ['NA', 'Undetermined', 'nan', 'None', '', '符合标准', '不符合标准']
                        else str(Decimal(str(float(x))).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                    )
            
            group_data = pd.concat([group_data, processed_group_result], ignore_index=True)
            all_results = pd.concat([all_results, group_data], ignore_index=True)
            
            # 添加空行（除了最后一组）
            if i < len(result) - 1:
                all_results = pd.concat([all_results, pd.DataFrame([['']], columns=['Sample Name'])], ignore_index=True)
    else:
        # 删除全是NA的列
        valid_cols = ['Sample Name']  # 始终保留Sample Name列
        for col in result.columns:
            if col != 'Sample Name':
                # 检查列中是否所有值都是NA或Undetermined
                values = result[col].astype(str)
                if not all(values.isin(['NA', 'Undetermined', 'nan', 'None', ''])):
                    valid_cols.append(col)
        result = result[valid_cols]
        
        # 添加数据（包括表头）
        result_data = pd.DataFrame(columns=result.columns)
        result_data.loc[0] = result.columns  # 添加表头
        
        # 处理数值数据
        processed_result = result.copy()
        for col in processed_result.columns:
            if col != 'Sample Name':
                processed_result[col] = processed_result[col].apply(
                    lambda x: str(x) if pd.isna(x) or not isinstance(x, (int, float)) or x in ['NA', 'Undetermined', 'nan', 'None', '', '符合标准', '不符合标准']
                    else str(Decimal(str(float(x))).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
                )
        
        result_data = pd.concat([result_data, processed_result], ignore_index=True)
        all_results = pd.concat([all_results, result_data], ignore_index=True)
        
        # 添加空行（除了最后一个模块）
        if not is_last_module:
            all_results = pd.concat([all_results, pd.DataFrame([['']], columns=['Sample Name'])], ignore_index=True)
    
    return all_results

def get_excel_formats(workbook):
    """获取Excel格式设置"""
    formats = {
        'title': workbook.add_format({
            'bold': True,
            'font_size': 14,
            'font_name': '微软雅黑',
            'align': 'center',
            'valign': 'vcenter',
            'bg_color': '#D9E1F2',
            'border': 1,
            'text_wrap': True
        }),
        'subtitle': workbook.add_format({
            'bold': True,
            'font_size': 12,
            'font_name': '微软雅黑',
            'align': 'center',
            'valign': 'vcenter',
            'bg_color': '#E2EFDA',
            'border': 1,
            'text_wrap': True
        }),
        'header': workbook.add_format({
            'bold': True,
            'font_size': 11,
            'font_name': '微软雅黑',
            'align': 'center',
            'valign': 'vcenter',
            'bg_color': '#F2F2F2',
            'border': 1,
            'text_wrap': True
        }),
        'data': workbook.add_format({
            'font_size': 10,
            'font_name': '微软雅黑',
            'align': 'center',
            'valign': 'vcenter',
            'border': 1,
            'text_wrap': True
        }),
        'empty': workbook.add_format({
            'top': 0,
            'bottom': 0,
            'left': 0,
            'right': 0
        })
    }
    return formats

def apply_excel_formats(worksheet, all_results, formats):
    """应用Excel格式"""
    for idx, row in all_results.iterrows():
        if row['Sample Name'].startswith('='):
            # 模块标题
            worksheet.set_row(idx, 30, formats['title'])
            # 合并单元格
            worksheet.merge_range(idx, 0, idx, len(all_results.columns)-1, row['Sample Name'], formats['title'])
        elif row['Sample Name'].startswith('-'):
            # 组标题
            worksheet.set_row(idx, 25, formats['subtitle'])
            # 合并单元格
            worksheet.merge_range(idx, 0, idx, len(all_results.columns)-1, row['Sample Name'], formats['subtitle'])
        elif row['Sample Name'] == '':
            # 空行
            worksheet.set_row(idx, 10, formats['empty'])
            # 清除空行的所有单元格边框
            for col in range(len(all_results.columns)):
                worksheet.write(idx, col, '', formats['empty'])
        else:
            # 检查是否是表头行
            is_header = False
            for col in all_results.columns:
                if col == row[col]:
                    is_header = True
                    break
            
            if is_header:
                # 表头行
                worksheet.set_row(idx, 20, formats['header'])
                # 写入表头
                for col_idx, col_name in enumerate(all_results.columns):
                    worksheet.write(idx, col_idx, col_name, formats['header'])
            else:
                # 数据行
                worksheet.set_row(idx, 20, formats['data'])
                # 写入数据
                for col_idx, col_name in enumerate(all_results.columns):
                    value = row[col_name]
                    # 处理NaN/INF值
                    if pd.isna(value) or (isinstance(value, float) and (np.isinf(value) or np.isnan(value))):
                        value = 'NA'
                    worksheet.write(idx, col_idx, value, formats['data'])

def export_comprehensive_report(analysis_results):
    """导出综合报告"""
    # 创建Excel文件
    excel_buffer = io.BytesIO()
    with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
        # 设置workbook选项
        workbook = writer.book
        workbook.nan_inf_to_errors = True
        
        # 创建一个空的DataFrame来存储所有结果
        all_results = pd.DataFrame()
        
        # 处理每个模块的数据
        modules = list(analysis_results.items())
        for i, (analysis_name, result) in enumerate(modules):
            if result is not None:
                module_results = process_module_data(
                    result, 
                    analysis_name, 
                    is_last_module=(i == len(modules) - 1)
                )
                all_results = pd.concat([all_results, module_results], ignore_index=True)
        
        # 将所有结果写入同一个sheet
        all_results.to_excel(writer, sheet_name='QPCR分析综合报告', index=False)
        
        # 获取工作表
        worksheet = writer.sheets['QPCR分析综合报告']
        
        # 获取格式设置
        formats = get_excel_formats(workbook)
        
        # 应用格式
        apply_excel_formats(worksheet, all_results, formats)
        
        # 调整列宽
        for idx, col in enumerate(all_results.columns):
            max_length = max(
                all_results[col].astype(str).apply(len).max(),
                len(str(col))
            )
            worksheet.set_column(idx, idx, max_length + 4)
    
    return excel_buffer

def show_comprehensive_export(analysis_results):
    """显示综合导出部分"""
    has_results = any(result is not None for result in analysis_results.values())
    if has_results:
        st.markdown("---")
        st.subheader("综合导出")
        
        # 导出综合报告
        excel_buffer = export_comprehensive_report(analysis_results)
        excel_buffer.seek(0)
        
        # 添加导出按钮
        st.download_button(
            label="导出综合报告(Excel)",
            data=excel_buffer,
            file_name=f"QPCR分析综合报告_{time.strftime('%Y%m%d_%H%M%S')}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
        
        # 添加说明
        st.info("综合报告包含以下内容：")
        for analysis_name, result in analysis_results.items():
            if result is not None:
                if analysis_name == '对照样品分析':
                    st.write(f"✓ {analysis_name}（{len(result)}组）")
                else:
                    st.write(f"✓ {analysis_name}")
            else:
                st.write(f"✗ {analysis_name}（未生成）")

def main():
    st.set_page_config(
        page_title="QPCR数据分析平台",
        page_icon="🧬",
        layout="wide"
    )
    
    # 简化标题
    st.title("QPCR数据分析平台")
    
    # 初始化分析结果存储
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = {
            '标准曲线分析': None,
            '质控样品分析': None,
            '对照样品分析': None,
            '待测样品分析': None,
            '汇报模块': None
        }
    
    # 添加参数管理部分
    st.sidebar.subheader("参数管理")
    
    # 创建参数存储字典
    if 'saved_params' not in st.session_state:
        st.session_state.saved_params = {}
    
    # 参数导出和导入功能
    with st.sidebar.expander("导出/导入参数"):
        # 导出参数
        if st.button("导出参数到文件"):
            current_params = {
                'saved_params': st.session_state.saved_params
            }
            
            # 创建参数文件
            param_buffer = io.BytesIO()
            pd.to_pickle(current_params, param_buffer)
            param_buffer.seek(0)
            
            # 提供下载
            st.download_button(
                label="下载参数文件",
                data=param_buffer,
                file_name=f"qpcr_params_{time.strftime('%Y%m%d_%H%M%S')}.pkl",
                mime="application/octet-stream"
            )
        
        # 导入参数
        uploaded_params = st.file_uploader("上传参数文件", type=['pkl'])
        if uploaded_params is not None:
            try:
                loaded_params = pd.read_pickle(uploaded_params)
                if 'saved_params' in loaded_params:
                    st.session_state.saved_params = loaded_params['saved_params']
                    st.success("参数加载成功！")
                else:
                    st.error("无效的参数文件格式")
            except Exception as e:
                st.error(f"参数加载失败: {str(e)}")
    
    # 参数保存功能
    with st.sidebar.expander("保存当前参数"):
        param_name = st.text_input("参数名称", key="save_param_name")
        if st.button("保存参数"):
            if param_name:
                current_params = {
                    'quantity_threshold': st.session_state.get('control_quantity_threshold', 0.5),
                    'theoretical_concentrations': st.session_state.get('theoretical_concentrations', {}),
                    'custom_thresholds': st.session_state.get('custom_thresholds', {}),
                    'target_genes': st.session_state.get('target_genes', {})
                }
                st.session_state.saved_params[param_name] = current_params
                st.success(f"参数 '{param_name}' 已保存")
            else:
                st.error("请输入参数名称")
    
    # 参数加载功能
    with st.sidebar.expander("加载已保存参数"):
        if st.session_state.saved_params:
            selected_param = st.selectbox(
                "选择要加载的参数",
                list(st.session_state.saved_params.keys())
            )
            if st.button("加载参数"):
                params = st.session_state.saved_params[selected_param]
                st.session_state.control_quantity_threshold = params['quantity_threshold']
                st.session_state.theoretical_concentrations = params['theoretical_concentrations']
                st.session_state.custom_thresholds = params['custom_thresholds']
                st.session_state.target_genes = params.get('target_genes', {})
                st.success(f"参数 '{selected_param}' 已加载")
        else:
            st.info("暂无保存的参数")
    
    # 参数删除功能
    with st.sidebar.expander("删除已保存参数"):
        if st.session_state.saved_params:
            param_to_delete = st.selectbox(
                "选择要删除的参数",
                list(st.session_state.saved_params.keys()),
                key="delete_param"
            )
            if st.button("删除参数"):
                del st.session_state.saved_params[param_to_delete]
                st.success(f"参数 '{param_to_delete}' 已删除")
        else:
            st.info("暂无保存的参数")
    
    # 1. 数据上传部分
    uploaded_file = st.file_uploader(
        "上传QPCR数据文件 (Excel格式)", 
        type=['xlsx', 'xls']
    )
    
    if uploaded_file is not None:
        try:
            # 从第47行开始读取数据
            df = pd.read_excel(uploaded_file, skiprows=46)
            df = df.dropna(how='all')  # 删除全为空的行
            
            # 2. 数据清洗
            target_dfs, cleaned_df = clean_qpcr_data(df)
            st.success("数据清洗完成")
            
            # 显示清洗后的数据
            st.write("数据清洗结果：")
            for target, target_df in target_dfs.items():
                with st.expander(f"查看 {target} 的清洗结果"):
                    st.dataframe(target_df)
            
            # 3. 分析模块
            tabs = st.tabs(["标准曲线分析", "质控样品分析", "对照样品分析", "待测样品分析", "汇报模块"])
            
            # 初始化目标基因参数存储
            if 'target_genes' not in st.session_state:
                st.session_state.target_genes = {}
            
            # 标准曲线分析
            with tabs[0]:
                st.subheader("标准曲线分析")
                
                col1, col2 = st.columns(2)
                with col1:
                    # 选择目标基因
                    target_names = list(target_dfs.keys())
                    selected_target = st.selectbox("选择目标基因", target_names)
                    
                    # 初始化当前目标基因的参数
                    if selected_target not in st.session_state.target_genes:
                        st.session_state.target_genes[selected_target] = {
                            'theoretical_concentrations': {},
                            'custom_thresholds': {}
                        }
                    
                    # 分别选择CV和RE的计算方式
                    cv_value_type = st.radio(
                        "选择CV计算方式",
                        ["CT值"],  # 只保留CT值选项
                        horizontal=True,
                        key="cv_type"
                    )
                    
                    re_value_type = st.radio(
                        "选择RE计算方式",
                        ["回算浓度"],
                        horizontal=True,
                        key="re_type"
                    )
                
                with col2:
                    # 样本选择
                    if selected_target:
                        target_df = target_dfs[selected_target]
                        all_samples = target_df['Sample Name'].unique().tolist()
                        
                        # 通过文本输入选择样本
                        st.write("选择要分析的样本：")
                        sample_input = st.text_area(
                            "输入要分析的样本",
                            height=150,
                            help="每行输入一个Sample Name"
                        )
                        
                        # 处理输入的样本
                        if sample_input:
                            # 分割输入文本并去除空行和空白字符
                            selected_samples = [s.strip() for s in sample_input.split('\n') if s.strip()]
                            
                            # 验证输入的Sample Name是否有效
                            invalid_samples = [s for s in selected_samples if s not in all_samples]
                            if invalid_samples:
                                st.error(f"以下Sample Name不存在：{', '.join(invalid_samples)}")
                                selected_samples = []
                        else:
                            selected_samples = all_samples  # 如果没有输入，使用所有样本
                
                # 添加每个样本的独立阈值设置
                with st.expander("设置样本独立阈值"):
                    custom_thresholds = st.session_state.target_genes[selected_target]['custom_thresholds']
                    samples_to_analyze = selected_samples if selected_samples else all_samples
                    
                    for sample in samples_to_analyze:
                        st.write(f"样本: {sample}")
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            cv = st.number_input(
                                "CV阈值 (%)",
                                value=custom_thresholds.get(sample, {}).get('cv', 5.0),
                                min_value=0.0,
                                step=0.1,
                                key=f"qc_cv_{selected_target}_{sample}_{int(time.time())}"
                            )
                        
                        with col2:
                            re = st.number_input(
                                "RE阈值 (%)",
                                value=custom_thresholds.get(sample, {}).get('re', 20.0),
                                min_value=0.0,
                                step=0.1,
                                key=f"qc_re_{selected_target}_{sample}_{int(time.time())}"
                            )
                        
                        if sample not in custom_thresholds:
                            custom_thresholds[sample] = {}
                        custom_thresholds[sample]['cv'] = cv
                        custom_thresholds[sample]['re'] = re
                
                # 设置理论浓度
                with st.expander("设置理论浓度"):
                    theoretical_concentrations = st.session_state.target_genes[selected_target]['theoretical_concentrations']
                    samples_to_analyze = selected_samples if selected_samples else all_samples
                    
                    for sample in samples_to_analyze:
                        theoretical_concentrations[sample] = st.number_input(
                            f"样本 {sample} 的理论浓度",
                            value=theoretical_concentrations.get(sample, 0.0),
                            min_value=0.0,
                            step=0.001,
                            format="%.3f",
                            key=f"qc_theo_{selected_target}_{sample}_{int(time.time())}"
                        )
                
                # 分析按钮
                if st.button("开始标准曲线分析", key=f"analyze_std_{selected_target}"):
                    if selected_target:
                        std_results, std_fig = analyze_standard_curve(
                            target_dfs[selected_target],
                            selected_target,
                            selected_samples,
                            cv_value_type,
                            re_value_type,
                            5.0,  # 使用固定值替代默认CV阈值
                            20.0,  # 使用固定值替代默认RE阈值
                            custom_thresholds=custom_thresholds,
                            theoretical_concentrations=theoretical_concentrations
                        )
                        
                        if std_results is not None:
                            st.write("标准曲线分析结果：")
                            st.dataframe(std_results)
                            
                            # 存储结果
                            st.session_state.analysis_results['标准曲线分析'] = std_results
                            
                            # 显示图表
                            if std_fig is not None:
                                st.plotly_chart(std_fig, key=f"std_curve_{selected_target}_{cv_value_type}_{re_value_type}_{int(time.time())}")
                            
                            # 添加导出功能
                            col1, col2 = st.columns(2)
                            with col1:
                                # 导出Excel
                                excel_buffer = io.BytesIO()
                                with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                                    std_results.to_excel(writer, sheet_name='标准曲线分析结果', index=False)
                                excel_buffer.seek(0)
                                st.download_button(
                                    label="导出Excel",
                                    data=excel_buffer,
                                    file_name=f"标准曲线分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.xlsx",
                                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                                )
                            
                            with col2:
                                # 导出CSV
                                csv_buffer = io.StringIO()
                                std_results.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
                                csv_buffer.seek(0)
                                st.download_button(
                                    label="导出CSV",
                                    data=csv_buffer.getvalue(),
                                    file_name=f"标准曲线分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                                    mime="text/csv"
                                )
                
                # 显示已保存的结果
                if st.session_state.analysis_results['标准曲线分析'] is not None:
                    st.write("已保存的标准曲线分析结果：")
                    st.dataframe(st.session_state.analysis_results['标准曲线分析'])
            
            # 质控样品分析
            with tabs[1]:
                st.subheader("质控样品分析")
                
                col1, col2 = st.columns(2)
                with col1:
                    # 选择目标基因
                    target_names = list(target_dfs.keys())
                    selected_target = st.selectbox("选择目标基因", target_names, key='qc_target')
                    
                    # 初始化当前目标基因的参数
                    if selected_target not in st.session_state.target_genes:
                        st.session_state.target_genes[selected_target] = {
                            'theoretical_concentrations': {},
                            'custom_thresholds': {}
                        }
                    
                    # 分别选择CV和RE的计算方式
                    cv_value_type = st.radio(
                        "选择CV计算方式",
                        ["CT值", "Quantity"],
                        horizontal=True,
                        key="qc_cv_type"
                    )
                    
                    re_value_type = st.radio(
                        "选择RE计算方式",
                        ["CT值", "Quantity"],
                        horizontal=True,
                        key="qc_re_type"
                    )
                
                with col2:
                    if selected_target:
                        target_df = target_dfs[selected_target]
                        all_samples = target_df['Sample Name'].unique().tolist()
                        
                        # 通过文本输入选择样本
                        st.write("选择要分析的样本：")
                        sample_input = st.text_area(
                            "输入要分析的样本",
                            height=150,
                            help="每行输入一个Sample Name",
                            key='qc_samples_input'
                        )
                        
                        # 处理输入的样本
                        if sample_input:
                            selected_samples = [s.strip() for s in sample_input.split('\n') if s.strip()]
                            invalid_samples = [s for s in selected_samples if s not in all_samples]
                            if invalid_samples:
                                st.error(f"以下Sample Name不存在：{', '.join(invalid_samples)}")
                                selected_samples = []
                        else:
                            selected_samples = all_samples
                
                # 添加每个样本的独立阈值设置
                with st.expander("设置样本独立阈值"):
                    custom_thresholds = st.session_state.target_genes[selected_target]['custom_thresholds']
                    samples_to_analyze = selected_samples if selected_samples else all_samples
                    
                    for sample in samples_to_analyze:
                        st.write(f"样本: {sample}")
                        col1, col2 = st.columns(2)
                        
                        # 生成随机数作为key的一部分
                        random_key = np.random.randint(1000000, 9999999)
                        
                        with col1:
                            cv = st.number_input(
                                "CV阈值 (%)",
                                value=custom_thresholds.get(sample, {}).get('cv', 5.0),
                                min_value=0.0,
                                step=0.1,
                                key=f"qc_cv_{selected_target}_{sample}_{random_key}"
                            )
                        
                        with col2:
                            re = st.number_input(
                                "RE阈值 (%)",
                                value=custom_thresholds.get(sample, {}).get('re', 20.0),
                                min_value=0.0,
                                step=0.1,
                                key=f"qc_re_{selected_target}_{sample}_{random_key}"
                            )
                        
                        if sample not in custom_thresholds:
                            custom_thresholds[sample] = {}
                        custom_thresholds[sample]['cv'] = cv
                        custom_thresholds[sample]['re'] = re
                
                # 设置理论浓度
                with st.expander("设置理论浓度"):
                    theoretical_concentrations = st.session_state.target_genes[selected_target]['theoretical_concentrations']
                    samples_to_analyze = selected_samples if selected_samples else all_samples
                    
                    for sample in samples_to_analyze:
                        # 生成随机数作为key的一部分
                        random_key = np.random.randint(1000000, 9999999)
                        
                        theoretical_concentrations[sample] = st.number_input(
                            f"样本 {sample} 的理论浓度",
                            value=theoretical_concentrations.get(sample, 0.0),
                            min_value=0.0,
                            step=0.001,
                            format="%.3f",
                            key=f"qc_theo_{selected_target}_{sample}_{random_key}"
                        )
                
                # 分析按钮
                if st.button("开始质控分析", key=f"analyze_qc_{selected_target}"):
                    if selected_target:
                        qc_results, qc_fig = analyze_qc_samples(
                            target_dfs[selected_target],
                            selected_target,
                            selected_samples,
                            cv_value_type,
                            re_value_type,
                            5.0,  # 使用固定值替代默认CV阈值
                            20.0,  # 使用固定值替代默认RE阈值
                            custom_thresholds=custom_thresholds,
                            theoretical_concentrations=theoretical_concentrations
                        )
                        
                        if qc_results is not None:
                            st.write("质控分析结果：")
                            st.dataframe(qc_results)
                            
                            # 存储结果
                            st.session_state.analysis_results['质控样品分析'] = qc_results
                            
                            # 添加导出功能
                            col1, col2 = st.columns(2)
                            with col1:
                                # 导出Excel
                                excel_buffer = io.BytesIO()
                                with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                                    qc_results.to_excel(writer, sheet_name='质控分析结果', index=False)
                                excel_buffer.seek(0)
                                st.download_button(
                                    label="导出Excel",
                                    data=excel_buffer,
                                    file_name=f"质控分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.xlsx",
                                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                                )
                            
                            with col2:
                                # 导出CSV
                                csv_buffer = io.StringIO()
                                qc_results.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
                                csv_buffer.seek(0)
                                st.download_button(
                                    label="导出CSV",
                                    data=csv_buffer.getvalue(),
                                    file_name=f"质控分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                                    mime="text/csv"
                                )
                
                # 显示已保存的结果
                if st.session_state.analysis_results['质控样品分析'] is not None:
                    st.write("已保存的质控分析结果：")
                    st.dataframe(st.session_state.analysis_results['质控样品分析'])
            
            # 对照样品分析
            with tabs[2]:
                st.subheader("对照样品分析")
                
                # 选择目标基因
                target_names = list(target_dfs.keys())
                selected_target = st.selectbox("选择目标基因", target_names, key='control_target')
                
                if selected_target:
                    target_df = target_dfs[selected_target]
                    all_samples = target_df['Sample Name'].unique().tolist()
                    
                    # 添加多组样品选择功能
                    st.write("选择要分析的样品组：")
                    num_groups = st.number_input(
                        "设置样品组数量",
                        min_value=1,
                        max_value=10,
                        value=1,
                        key='control_num_groups'
                    )
                    
                    # 为每组样品创建选择区域和标准设置
                    selected_samples_groups = []
                    group_settings = []
                    
                    for i in range(num_groups):
                        st.markdown(f"### 第 {i+1} 组设置")
                        
                        # 样品选择
                        st.write("选择要分析的样本：")
                        sample_input = st.text_area(
                            f"输入第 {i+1} 组要分析的样本",
                            height=100,
                            help="每行输入一个Sample Name，可以输入不存在的Sample Name用于参数设置",
                            key=f'control_samples_input_{i}'
                        )
                        
                        # 处理输入的样本
                        samples = []
                        if sample_input:
                            samples = [s.strip() for s in sample_input.split('\n') if s.strip()]
                            # 检查是否存在无效样本，但只显示警告
                            invalid_samples = [s for s in samples if s not in all_samples]
                            if invalid_samples:
                                st.warning(f"注意：以下Sample Name不存在：{', '.join(invalid_samples)}")
                        
                        if samples:
                            selected_samples_groups.append(samples)
                            
                            # 为每组设置独立的标准
                            col1, col2 = st.columns(2)
                            with col1:
                                value_type = st.radio(
                                    "选择分析值类型",
                                    ["CT值", "Quantity"],
                                    horizontal=True,
                                    key=f"control_value_type_{i}"
                                )
                                
                                threshold_type = st.radio(
                                    "选择标准类型",
                                    ["小于", "小于等于", "大于", "大于等于", "范围"],
                                    horizontal=True,
                                    key=f"control_threshold_type_{i}"
                                )
                            
                            with col2:
                                if threshold_type == "范围":
                                    min_threshold = st.number_input(
                                        "最小值",
                                        value=st.session_state.get(f'control_min_threshold_{i}', 0.0),
                                        min_value=0.0,
                                        step=0.1,
                                        format="%.3f",
                                        key=f"control_min_threshold_{i}"
                                    )
                                    max_threshold = st.number_input(
                                        "最大值",
                                        value=st.session_state.get(f'control_max_threshold_{i}', 1.0),
                                        min_value=0.0,
                                        step=0.1,
                                        format="%.3f",
                                        key=f"control_max_threshold_{i}"
                                    )
                                else:
                                    threshold = st.number_input(
                                        "阈值",
                                        value=st.session_state.get(f'control_threshold_{i}', 0.5),
                                        min_value=0.0,
                                        step=0.1,
                                        format="%.3f",
                                        key=f"control_threshold_{i}"
                                    )
                            
                            # 存储该组的设置
                            group_settings.append({
                                'value_type': value_type,
                                'threshold_type': threshold_type,
                                'threshold': threshold if threshold_type != "范围" else None,
                                'min_threshold': min_threshold if threshold_type == "范围" else None,
                                'max_threshold': max_threshold if threshold_type == "范围" else None
                            })
                    
                    # 分析数据
                    if st.button("开始分析", key=f"control_analyze_{selected_target}"):
                        if selected_samples_groups and group_settings:
                            all_results = []
                            for i, (samples, settings) in enumerate(zip(selected_samples_groups, group_settings)):
                                st.write(f"第 {i+1} 组分析结果：")
                                # 过滤出存在的样本进行分析
                                valid_samples = [s for s in samples if s in all_samples]
                                if valid_samples:
                                    results = analyze_control_samples(
                                        target_df,
                                        valid_samples,
                                        settings['value_type'],
                                        settings['threshold_type'],
                                        settings['threshold'],
                                        settings['min_threshold'],
                                        settings['max_threshold']
                                    )
                                    
                                    if results is not None:
                                        st.dataframe(results)
                                        all_results.append(results)
                                else:
                                    st.warning(f"第 {i+1} 组没有有效的样本可供分析。")
                            
                            # 存储结果
                            st.session_state.analysis_results['对照样品分析'] = all_results
                            
                            # 添加导出功能
                            if all_results:
                                col1, col2 = st.columns(2)
                                with col1:
                                    # 导出Excel
                                    excel_buffer = io.BytesIO()
                                    with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                                        for i, result in enumerate(all_results):
                                            result.to_excel(writer, sheet_name=f'对照分析结果_组{i+1}', index=False)
                                    excel_buffer.seek(0)
                                    st.download_button(
                                        label="导出Excel",
                                        data=excel_buffer,
                                        file_name=f"对照分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.xlsx",
                                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                                    )
                                
                                with col2:
                                    # 导出CSV
                                    csv_buffer = io.StringIO()
                                    for i, result in enumerate(all_results):
                                        csv_buffer.write(f"第{i+1}组分析结果\n")
                                        result.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
                                        csv_buffer.write("\n")
                                    csv_buffer.seek(0)
                                    st.download_button(
                                        label="导出CSV",
                                        data=csv_buffer.getvalue(),
                                        file_name=f"对照分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                                        mime="text/csv"
                                    )
                        else:
                            st.error("请至少选择一组样品进行分析。")
                
                # 显示已保存的结果
                if st.session_state.analysis_results['对照样品分析'] is not None:
                    st.write("已保存的对照分析结果：")
                    for i, result in enumerate(st.session_state.analysis_results['对照样品分析']):
                        st.write(f"第 {i+1} 组结果：")
                        st.dataframe(result)
            
            # 待测样品分析
            with tabs[3]:
                st.subheader("待测样品分析")
                
                # 选择目标基因
                target_names = list(target_dfs.keys())
                selected_target = st.selectbox("选择目标基因", target_names, key='test_target')
                
                if selected_target:
                    target_df = target_dfs[selected_target]
                    all_samples = target_df['Sample Name'].unique().tolist()
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        # 选择分析值类型
                        value_type = st.radio(
                            "选择分析值类型",
                            ["CT值", "Quantity", "CT值和Quantity"],
                            horizontal=True,
                            key="test_value_type"
                        )
                        
                        # 选择是否进行CV分析
                        do_cv_analysis = st.checkbox(
                            "进行CV分析",
                            value=True,
                            key="test_cv_check"
                        )
                        
                        # 选择是否进行阈值分析
                        do_threshold_analysis = st.checkbox(
                            "进行阈值分析",
                            value=False,
                            key="test_threshold_check"
                        )
                    
                    with col2:
                        # 通过文本输入选择样本
                        st.write("选择要分析的样本：")
                        sample_input = st.text_area(
                            "输入要分析的样本",
                            height=150,
                            help="每行输入一个Sample Name",
                            key='test_samples_input'
                        )
                        
                        # 处理输入的样本
                        if sample_input:
                            selected_samples = [s.strip() for s in sample_input.split('\n') if s.strip()]
                            invalid_samples = [s for s in selected_samples if s not in all_samples]
                            if invalid_samples:
                                st.error(f"以下Sample Name不存在：{', '.join(invalid_samples)}")
                                selected_samples = []
                        else:
                            selected_samples = all_samples
                    
                    # 如果选择进行CV分析，显示CV设置
                    cv_threshold = None
                    if do_cv_analysis:
                        with st.expander("设置CV接受标准"):
                            cv_threshold = st.number_input(
                                "CV接受标准 (%)",
                                value=5.0,
                                min_value=0.0,
                                step=0.1,
                                format="%.1f",
                                key="test_cv_threshold"
                            )
                    
                    # 如果选择进行阈值分析，显示阈值设置
                    threshold_type = None
                    threshold = None
                    min_threshold = None
                    max_threshold = None
                    
                    if do_threshold_analysis:
                        with st.expander("设置阈值标准"):
                            threshold_type = st.radio(
                                "选择标准类型",
                                ["小于", "小于等于", "大于", "大于等于", "范围"],
                                horizontal=True,
                                key="test_threshold_type"
                            )
                            
                            if threshold_type == "范围":
                                col1, col2 = st.columns(2)
                                with col1:
                                    min_threshold = st.number_input(
                                        "最小值",
                                        value=0.0,
                                        min_value=0.0,
                                        step=0.1,
                                        format="%.3f",
                                        key="test_min_threshold"
                                    )
                                with col2:
                                    max_threshold = st.number_input(
                                        "最大值",
                                        value=1.0,
                                        min_value=0.0,
                                        step=0.1,
                                        format="%.3f",
                                        key="test_max_threshold"
                                    )
                            else:
                                threshold = st.number_input(
                                    "阈值",
                                    value=0.5,
                                    min_value=0.0,
                                    step=0.1,
                                    format="%.3f",
                                    key="test_threshold"
                                )
                    
                    # 分析按钮
                    if st.button("开始分析", key=f"analyze_test_{selected_target}"):
                        if selected_target:
                            test_results, test_fig = analyze_test_samples(
                                target_dfs[selected_target],
                                selected_target,
                                selected_samples,
                                value_type,
                                threshold_type if do_threshold_analysis else None,
                                threshold if do_threshold_analysis else None,
                                min_threshold if do_threshold_analysis else None,
                                max_threshold if do_threshold_analysis else None,
                                cv_threshold if do_cv_analysis else None
                            )
                            
                            if test_results is not None:
                                st.write("待测样品分析结果：")
                                st.dataframe(test_results)
                                
                                # 存储结果
                                st.session_state.analysis_results['待测样品分析'] = test_results
                                
                                # 添加导出功能
                                col1, col2 = st.columns(2)
                                with col1:
                                    # 导出Excel
                                    excel_buffer = io.BytesIO()
                                    with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                                        test_results.to_excel(writer, sheet_name='待测样品分析结果', index=False)
                                    excel_buffer.seek(0)
                                    st.download_button(
                                        label="导出Excel",
                                        data=excel_buffer,
                                        file_name=f"待测样品分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.xlsx",
                                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                                    )
                                
                                with col2:
                                    # 导出CSV
                                    csv_buffer = io.StringIO()
                                    test_results.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
                                    csv_buffer.seek(0)
                                    st.download_button(
                                        label="导出CSV",
                                        data=csv_buffer.getvalue(),
                                        file_name=f"待测样品分析结果_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                                        mime="text/csv"
                                    )
                    
                    # 显示已保存的结果
                    if st.session_state.analysis_results['待测样品分析'] is not None:
                        st.write("已保存的待测样品分析结果：")
                        st.dataframe(st.session_state.analysis_results['待测样品分析'])
            
            # 汇报模块
            with tabs[4]:
                st.subheader("汇报模块")
                
                # 选择目标基因
                target_names = list(target_dfs.keys())
                selected_target = st.selectbox("选择目标基因", target_names, key='report_target')
                
                if selected_target:
                    target_df = target_dfs[selected_target]
                    all_samples = target_df['Sample Name'].unique().tolist()
                    
                    # 选择汇报类型
                    report_type = st.radio(
                        "选择汇报类型",
                        ["VCN", "RCL", "DNA中拷贝数"],
                        horizontal=True,
                        key="report_type"
                    )
                    
                    # 通过文本输入选择样本
                    st.write("选择要分析的样本：")
                    sample_input = st.text_area(
                        "输入要分析的样本",
                        height=150,
                        help="每行输入一个Sample Name",
                        key='report_samples_input'
                    )
                    
                    # 处理输入的样本
                    if sample_input:
                        selected_samples = [s.strip() for s in sample_input.split('\n') if s.strip()]
                        invalid_samples = [s for s in selected_samples if s not in all_samples]
                        if invalid_samples:
                            st.error(f"以下Sample Name不存在：{', '.join(invalid_samples)}")
                            selected_samples = []
                    else:
                        selected_samples = all_samples
                    
                    # 如果是RCL汇报，添加阴阳性判定规则设置
                    if report_type == "RCL":
                        with st.expander("设置阴阳性判定规则"):
                            # 设置CT阈值
                            threshold_value = st.number_input(
                                "CT阈值",
                                value=40.0,
                                min_value=0.0,
                                step=0.1,
                                format="%.3f",
                                key="determination_threshold"
                            )
                            st.info("判定规则：CT值 > 阈值为阴性，CT值 ≤ 阈值为阳性")
                    
                    if st.button("生成汇报", key="generate_report"):
                        report_data = generate_report(
                            target_dfs,
                            selected_target,
                            selected_samples if selected_samples else None,
                            report_type,
                            "CT值",  # 固定使用CT值
                            "小于等于",  # 固定使用小于等于
                            threshold_value if report_type == "RCL" else None,
                            None,
                            None
                        )
                        
                        if report_data is not None and not report_data.empty:
                            st.write(f"{report_type}汇报结果：")
                            st.dataframe(report_data)
                            
                            # 存储结果
                            st.session_state.analysis_results['汇报模块'] = report_data
                            
                            # 添加导出功能
                            col1, col2 = st.columns(2)
                            with col1:
                                # 导出Excel
                                excel_buffer = io.BytesIO()
                                with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                                    report_data.to_excel(writer, sheet_name=f'{report_type}汇报数据', index=False)
                                excel_buffer.seek(0)
                                st.download_button(
                                    label="导出Excel",
                                    data=excel_buffer,
                                    file_name=f"{report_type}汇报数据_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.xlsx",
                                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                                )
                            
                            with col2:
                                # 导出CSV
                                csv_buffer = io.StringIO()
                                report_data.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
                                csv_buffer.seek(0)
                                st.download_button(
                                    label="导出CSV",
                                    data=csv_buffer.getvalue(),
                                    file_name=f"{report_type}汇报数据_{selected_target}_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                                    mime="text/csv"
                                )
                    
                    # 显示已保存的结果
                    if st.session_state.analysis_results['汇报模块'] is not None:
                        st.write("已保存的汇报结果：")
                        st.dataframe(st.session_state.analysis_results['汇报模块'])
            
            # 添加综合导出功能
            show_comprehensive_export(st.session_state.analysis_results)
        
        except Exception as e:
            st.error(f"数据处理出错: {str(e)}")

if __name__ == "__main__":
    main() 
