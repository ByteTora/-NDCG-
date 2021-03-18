import math
import numpy as np
import pandas as pd
import itertools

# 一. 数据准备
df_data = pd.read_excel('./original_table.xlsx')  # 导入数据表
name = df_data['小分子名称']  # 记录表头
name2IC50 = dict(zip(name, df_data['IC50值']))  # 将“小分子名称“和”IC50值“对应存入一个新字典
cycle_index = len(name2IC50)  # 计算小分子个数，为取数做准备


# 二. 编写计算IDCG和DCG值的函数

## 1. 计算IC50值的DCG或IDCG
def calculate_IDCG_or_DCG(rels):
    caculate_output = 0
    rels = [-math.log10(rel) for rel in rels]
    for i in range(1, len(rels)):
        caculate_output += rels[i] / math.log2(i + 1)
    caculate_output = caculate_output + rels[0]
    return caculate_output


## 2. 计算IC50的nDCG
def calculate_lis_nDCG(all_combinations_lis):
    lis_ndcg = []
    for single_combination in all_combinations_lis:
        lis_IC50 = []
        lis_sorted_IC50 = []
        for i in single_combination:
            lis_IC50.append(name2IC50[i])  # 取小分子名称对应的IC50值
            lis_sorted_IC50.append(name2IC50[i])
        lis_sorted_IC50 = sorted(lis_sorted_IC50)  # 将未取-lg的IC50值升序排列
        # 调用函数计算IC50的DCG的值
        dcg = calculate_IDCG_or_DCG(lis_IC50)
        # 调用函数计算IC50的IDCG的值
        idcg = calculate_IDCG_or_DCG(lis_sorted_IC50)
        ndcg = dcg / idcg
        lis_ndcg.append(ndcg)
    return lis_ndcg


## 3. 取数函数
def take_number_from_table_com(num):
    # 此方法的作用是:从name个小分子中取num个小分子的所有组合情况, 即C(name,num)
    iter = itertools.combinations(name, num)
    return iter


# 三. 主函数
if __name__ == '__main__':
    max_ndcg_number = 0
    for j in range(1, cycle_index):  # 从小分子中每次取j个小分子
        # 将从j个小分子中取i个小分子的所有组合情况存入一个列表
        single_all_combinations_lis = list(take_number_from_table_com(j + 1))
        # 计算j个小分子中取i个小分子的所有组合情况的nDCG的索引，并存入一个列表
        single_all_ndcg_lis = calculate_lis_nDCG(single_all_combinations_lis)
        # 调用numpy中的max方法，计算j个小分子中取i个小分子的所有组合情况的最大nDCG
        single_max_ncdg_number = np.max(single_all_ndcg_lis)
        if max_ndcg_number < single_max_ncdg_number:  # 判断此次（j小分子中取i个小分子所有情况的最大nDCG）是否为所有小分子的最大的nDCG
            max_ndcg_number = single_max_ncdg_number
            # 调用numpy的argmax方法，求出所有情况的最大nDCG对应的索引（即小分子序列）
            max_combinations_lis = single_all_combinations_lis[np.argmax(single_all_ndcg_lis)]
    print(f"所有小分子的所有组合情况中使得nDCG最大的值的序列为{max_combinations_lis}，其值为{max_ndcg_number}")
