# import numpy as np
#
#
# k=3
# a=np.array([7,5,1,4,9,6,3,4,8,6,5,8])
# b = a.argsort()[-k:]
# print(b)
# print(list(a[list(b)]))
import math
import statistics

import numpy as np

# a=float('3.5')
# b=np.array(a).astype(np.float32).mean()
# c=[b,b]
# d=statistics.mean(c)
# print(d)
# import random
# real_time=[]
# for i in range( 1):
#     rate = random.uniform(-0.2, 0.2)
#     real_time.append(20000 * (1 + rate))
# time_serial = sum(real_time)
# time_parallel = max(real_time)
# print(time_parallel)
# print(time_serial)
# result_fragment=['C=O','c1ccccc1']
# scores=['1','2']


# from rdkit import Chem
# from rdkit.Chem import Draw
# png_path =  'test.png'
# img = Draw.MolsToGridImage(mols=[Chem.MolFromSmiles(x) for x in result_fragment], legends=scores)
# img.save(png_path)

# a=[1,2]
# b=':'.join(a)
# print(b)
# import argparse
# parser = argparse.ArgumentParser(description='BO Search with Chemprop')
# parser.add_argument('--test',action='store_true', help='config of SMINA')
# args = parser.parse_args()
#
# print(args.test)

# def adaptive_function(iter,num_sampling,alpha,b):
#     return num_sampling*math.exp(-1*alpha*(iter-1))+b*math.log(iter)
#
# for i in range(1,8):
#     m=adaptive_function(i,8000,1,100)
#     print(int(m))
# new_mol='这是一个新的有效分子结构，类似于C1(CC2=C(C1)C(=O)NC2)Nc3ccc5c(c3)CC[C@@H]5：```C1(CC2=C(C1)C(=O)NC2)Nc3cc4c(c3)CC[C@@H]4```'
# while True:
#     filter = r'\'(.*?)\''
#     import re
#
#     new_mol_filted = re.findall(filter, new_mol)
#     if len(new_mol_filted) != 0:
#         new_mol = new_mol_filted[0]
#     else:
#         break
# print(new_mol)
# import numpy as np
# population=np.array([1,2,3,4,5,6,7,8,9,10])
# weights=np.array([10,10,10,10,10,10.1,10.1,10.2,10.2,10.2])
# for i in range(10):
#     sample = np.random.choice(population, size=3, replace=False, p=weights / sum(weights))
#     print(sample)
# a=[1,2,3,4,5,6,7,-1]
# print(a[-1])
# a.sort()
# print(a)



# ###begin
# import pandas as pd
# import rdkit
# from rdkit import Chem
# from tools.docking import calc_affinity
# result_file='/home/th2024/lmlf/LMLF-main/TSNE_8UOB_ours/result.csv'
#
# df=pd.read_csv(result_file)
# smiles=list(df.iloc[:,0])
# scores=list(df.iloc[:,1])
# print(scores[740])
# print(smiles[740])
# smiles1=[]
# scores1=[]
# for i in range(len(smiles)):
#     if scores[i]>9:
#         smiles1.append(smiles[i])
#         scores1.append(scores[i])
# print(len((smiles1)))
# for i in range(len((smiles1))):
#     if smiles1[i]=='CC1=C(NC(=O)C2=CC=CC(C)=C2C2=NC(N)=CC=C2)C=CC=C1N':
#         print(i)
# print(smiles1[24])
# print(scores1[24])
# from rdkit.Chem import RDConfig
# import sys
# import os
# sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# import sascorer
# m = Chem.MolFromSmiles(smiles1[24])
# sa_score = sascorer.calculateScore(m)
# print(sa_score)
# SA2=[]
# smiles2=[]
# scores2=[]
# for i in range(len(smiles1)):
#     m = Chem.MolFromSmiles(smiles1[i])
#     sa_score = sascorer.calculateScore(m)
#     if sa_score<4:
#         smiles2.append(smiles1[i])
#         scores2.append(scores1[i])
#         SA2.append(sa_score)
# print(len(smiles2))
#
# for i in range(len(smiles2)):
#     name='result_'+str(i)+'_'
#     smile=smiles2[i]
#     calc_affinity(sml=smile,dir_out='8UOB_result_filted',prefix=name)
# ###end


# xvg_to_csv.py
# import os
# def xvg_to_csv(input_file, output_file):
#    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
#        for line in infile:
#            if not line.startswith(('#', '@')): # 跳过注释和元数据行
#                outfile.write(line.replace(' ', ',')) # 替换空格为逗号
#
#
# path='/home/th2024/8UOB_gromacs/'
# name='hbnum'
# input_name=name+'.xvg'
# out_name=name+'.csv'
# input_path=os.path.join(path,input_name)
# output_path=os.path.join(path,out_name)
# xvg_to_csv(input_path,output_path)

# !/usr/bin/env python3
"""
将GROMACS的.xvg文件转换为.csv文件
支持自动检测注释行、多列数据、以及不同的数据格式
"""

import pandas as pd
import numpy as np
import re
import sys
import os
from pathlib import Path


def parse_xvg_file(xvg_file, delimiter=None):
    """
    解析.xvg文件，提取数据和元数据

    参数:
    xvg_file: .xvg文件路径
    delimiter: 列分隔符（如为None则自动检测）

    返回:
    dict: 包含数据和元数据的字典
    """
    data = {
        'metadata': {},
        'comments': [],
        'column_names': [],
        'data': []
    }

    # 自动检测常见分隔符
    possible_delimiters = [None, ' ', '\t', ',', ';']

    with open(xvg_file, 'r') as f:
        lines = f.readlines()

    # 解析文件
    for i, line in enumerate(lines):
        line = line.strip()

        # 空行
        if not line:
            continue

        # 注释行（以@或#开头）
        if line.startswith('#') or line.startswith('@'):
            data['comments'].append(line)

            # 尝试提取列名（从@ s? legend "Column Name"格式）
            if 'legend' in line:
                match = re.search(r'legend\s+"([^"]+)"', line)
                if match:
                    data['column_names'].append(match.group(1))

            # 提取其他元数据
            if line.startswith('@'):
                # 例如: @ title "RMSD"
                match = re.search(r'@\s+(\w+)\s+"([^"]+)"', line)
                if match:
                    key, value = match.groups()
                    data['metadata'][key] = value
            continue

        # 数据行（第一个非注释行通常是列数信息）
        if i == len([l for l in lines if l.strip() and not (l.startswith('#') or l.startswith('@'))]):
            # 这可能是标题行或列数声明
            if 'xaxis' in line.lower() or 'yaxis' in line.lower():
                continue

        # 真正的数据行
        if delimiter is None:
            # 尝试自动检测分隔符
            for delim in ['\t', ' ', ',', ';']:
                parts = line.split(delim)
                if len(parts) > 1 and all(self._is_number(p) for p in parts[:min(2, len(parts))]):
                    delimiter = delim
                    break

        # 分割数据
        if delimiter:
            parts = line.split(delimiter)
        else:
            # 默认按空白字符分割（兼容多个空格或制表符）
            parts = re.split(r'\s+', line)

        # 清理和转换数据
        clean_parts = []
        for part in parts:
            if part:  # 跳过空字符串
                try:
                    # 尝试转换为浮点数
                    clean_parts.append(float(part))
                except ValueError:
                    # 如果不能转换，保持原样
                    clean_parts.append(part)

        if clean_parts:
            data['data'].append(clean_parts)

    # 确保所有数据行长度一致
    if data['data']:
        max_len = max(len(row) for row in data['data'])
        for i in range(len(data['data'])):
            if len(data['data'][i]) < max_len:
                data['data'][i].extend([np.nan] * (max_len - len(data['data'][i])))

    return data


def _is_number(s):
    """检查字符串是否可以转换为数字"""
    try:
        float(s)
        return True
    except ValueError:
        return False


def convert_xvg_to_csv(xvg_file, csv_file=None, include_comments=False,
                       custom_headers=None, delimiter=','):
    """
    主函数：将.xvg文件转换为.csv文件

    参数:
    xvg_file: 输入.xvg文件路径
    csv_file: 输出.csv文件路径（如为None，则自动生成）
    include_comments: 是否在CSV中包含注释作为前几行
    custom_headers: 自定义列名列表
    delimiter: CSV文件使用的分隔符（默认逗号）
    """

    # 输入文件检查
    if not os.path.exists(xvg_file):
        print(f"错误: 文件不存在 - {xvg_file}")
        return False

    # 输出文件名
    if csv_file is None:
        csv_file = Path(xvg_file).with_suffix('.csv')
    else:
        csv_file = Path(csv_file)

    # 解析.xvg文件
    print(f"正在解析文件: {xvg_file}")
    parsed_data = parse_xvg_file(xvg_file)

    if not parsed_data['data']:
        print("警告: 未找到数据行")
        return False

    # 创建DataFrame
    data_array = np.array(parsed_data['data'])
    df = pd.DataFrame(data_array)

    # 设置列名
    if custom_headers and len(custom_headers) == df.shape[1]:
        df.columns = custom_headers
    elif parsed_data['column_names'] and len(parsed_data['column_names']) == df.shape[1]:
        df.columns = parsed_data['column_names']
    else:
        # 默认列名
        df.columns = [f'Column_{i + 1}' for i in range(df.shape[1])]
        # 通常第一列是时间，第二列是主要数据
        if df.shape[1] >= 2:
            df.columns = ['Time'] + [f'Data_{i}' for i in range(df.shape[1] - 1)]

    # 保存为CSV
    try:
        if include_comments and parsed_data['comments']:
            # 先写入注释
            with open(csv_file, 'w') as f:
                for comment in parsed_data['comments']:
                    f.write(f"# {comment}\n")
                # 再写入数据
                df.to_csv(f, index=False, sep=delimiter)
            print(f"已保存带注释的CSV文件: {csv_file}")
        else:
            df.to_csv(csv_file, index=False, sep=delimiter)
            print(f"已保存CSV文件: {csv_file}")

        # 打印摘要信息
        print(f"数据摘要:")
        print(f"  数据行数: {len(df)}")
        print(f"  数据列数: {len(df.columns)}")
        print(f"  列名: {list(df.columns)}")
        print(f"  数据类型:\n{df.dtypes.to_string()}")

        return True

    except Exception as e:
        print(f"保存CSV时出错: {e}")
        return False


def batch_convert_xvg(input_dir, output_dir=None, pattern="*.xvg"):
    """
    批量转换目录中的所有.xvg文件

    参数:
    input_dir: 输入目录
    output_dir: 输出目录（如为None，则使用输入目录）
    pattern: 文件匹配模式
    """
    input_path = Path(input_dir)

    if output_dir is None:
        output_dir = input_dir
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # 查找所有.xvg文件
    xvg_files = list(input_path.glob(pattern))

    if not xvg_files:
        print(f"在 {input_dir} 中未找到 {pattern} 文件")
        return

    print(f"找到 {len(xvg_files)} 个.xvg文件，开始批量转换...")

    success_count = 0
    for xvg_file in xvg_files:
        csv_file = output_path / (xvg_file.stem + '.csv')
        print(f"\n处理文件: {xvg_file.name}")

        if convert_xvg_to_csv(xvg_file, csv_file):
            success_count += 1

    print(f"\n批量转换完成!")
    print(f"成功: {success_count}/{len(xvg_files)}")
    print(f"输出目录: {output_path}")


def main():
    """命令行入口点"""
    import argparse

    parser = argparse.ArgumentParser(
        description='将GROMACS的.xvg文件转换为.csv文件',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 转换单个文件
  python xvg_to_csv.py simulation.xvg

  # 转换单个文件并指定输出名
  python xvg_to_csv.py input.xvg -o output.csv

  # 批量转换目录中的所有.xvg文件
  python xvg_to_csv.py -d ./results/

  # 使用自定义列名
  python xvg_to_csv.py rmsd.xvg -c "Time" "RMSD (nm)"

  # 使用分号作为CSV分隔符
  python xvg_to_csv.py energy.xvg --csv-delimiter ';'
        """
    )

    parser.add_argument('--input', default='/home/th2024/8UOB_gromacs/rmsd.xvg', help='输入.xvg文件路径')
    parser.add_argument('-o', '--output',default='/home/th2024/8UOB_gromacs/rmsd.xvg', help='输出.csv文件路径')
    parser.add_argument('-d', '--directory', help='批量转换目录')
    parser.add_argument('-c', '--columns', nargs='+', help='自定义列名')
    parser.add_argument('--include-comments', action='store_true',
                        help='在CSV中包含注释行')
    parser.add_argument('--csv-delimiter', default=',',
                        help='CSV文件分隔符 (默认: ",")')
    parser.add_argument('--pattern', default='*.xvg',
                        help='批量转换时的文件匹配模式 (默认: "*.xvg")')

    args = parser.parse_args()

    if args.directory:
        # 批量转换模式
        batch_convert_xvg(args.directory, args.output, args.pattern)
    elif args.input:
        # 单个文件转换模式
        convert_xvg_to_csv(
            args.input,
            args.output,
            args.include_comments,
            args.columns,
            args.csv_delimiter
        )
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()






