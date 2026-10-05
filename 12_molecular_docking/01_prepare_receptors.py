# -*- coding: utf-8 -*-
"""
---
title: "分子对接：受体制备（受体蛋白预处理）"
description: "从PDB数据库下载铁死亡关键蛋白结构，使用AutoDockTools进行加氢、计算电荷、添加极性氢等预处理步骤，准备用于分子对接的受体文件"
input: "PDB IDs列表（GPX4: 2OBI, Nrf2: 4IQK, SLC7A11: 6EWV, ACSL4: 6VJ2, HMOX1: 1N45）"
output: "results/15_docking/receptors/*.pdbqt （预处理后的受体文件）"
dependencies: "biopython, requests, pymol（可选，用于可视化）"
---

用法：
    python 01_prepare_receptors.py

依赖：
    - AutoDock Vina / AutoDockTools (MGLTools)
    - Biopython (可选，用于PDB文件解析)
    - 本脚本假设系统已安装MGLTools的prepare_receptor4.py脚本
"""

import os
import sys
import subprocess
import requests
from pathlib import Path


# ===================== 配置 =====================

# 项目根目录（从脚本位置向上一级）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 输出目录
RECEPTOR_DIR = PROJECT_ROOT / "results" / "15_docking" / "receptors"
PDB_RAW_DIR = PROJECT_ROOT / "results" / "15_docking" / "pdb_raw"

# 铁死亡关键蛋白 - PDB ID 映射
RECEPTORS = {
    "GPX4": "2OBI",      # Glutathione peroxidase 4
    "NFE2L2": "4IQK",    # Nrf2 (NFE2L2) - Kelch-like ECH-associated protein 1 complex
    "SLC7A11": "6EWV",   # xCT / System XC- light chain
    "ACSL4": "6VJ2",     # Acyl-CoA synthetase long-chain family member 4
    "HMOX1": "1N45",     # Heme oxygenase 1
    "PTGS2": "5IKQ",     # Cyclooxygenase-2 (COX-2)
    "TFRC": "3KAS",      # Transferrin receptor
    "HSPB1": "3Q9P"      # Heat shock protein beta-1
}

# MGLTools路径（根据实际安装位置修改）
MGLTOOLS_PATH = os.environ.get("MGLTOOLS_PATH", "/usr/local/MGLToolsPckgs")
PREPARE_RECEPTOR_SCRIPT = os.path.join(
    MGLTOOLS_PATH, "AutoDockTools/Utilities24/prepare_receptor4.py"
)


def download_pdb(pdb_id, output_dir):
    """从RCSB PDB下载PDB文件"""
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_file = output_dir / f"{pdb_id}.pdb"
    
    if output_file.exists():
        print(f"  [skip] {pdb_id}.pdb 已存在")
        return output_file
    
    print(f"  下载 {pdb_id} ...", end=" ")
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        output_file.write_text(response.text)
        print("完成")
        return output_file
    except Exception as e:
        print(f"失败: {e}")
        return None


def prepare_receptor(pdb_file, output_dir, protein_name):
    """
    使用AutoDockTools预处理受体：
    - 去水
    - 加极性氢
    - 加Gasteiger电荷
    - 保存为PDBQT格式
    """
    pdbqt_file = output_dir / f"{protein_name}.pdbqt"
    
    if pdbqt_file.exists():
        print(f"  [skip] {protein_name}.pdbqt 已存在")
        return pdbqt_file
    
    # 检查prepare_receptor4.py是否存在
    if not os.path.exists(PREPARE_RECEPTOR_SCRIPT):
        print(f"  警告: 未找到 {PREPARE_RECEPTOR_SCRIPT}")
        print(f"  跳过自动预处理，请手动使用AutoDockTools处理 {pdb_file}")
        print(f"  输出文件应为: {pdbqt_file}")
        return None
    
    cmd = [
        "pythonsh", PREPARE_RECEPTOR_SCRIPT,
        "-r", str(pdb_file),
        "-o", str(pdbqt_file),
        "-A", "hydrogens",  # 添加氢
        "-U", "waters",      # 去水
        "-v"
    ]
    
    print(f"  预处理 {protein_name} ...", end=" ")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if pdbqt_file.exists():
            print("完成")
            return pdbqt_file
        else:
            print(f"失败")
            print(f"  stdout: {result.stdout[-200:]}")
            print(f"  stderr: {result.stderr[-200:]}")
            return None
    except Exception as e:
        print(f"错误: {e}")
        return None


def main():
    """主函数：下载并预处理所有受体蛋白"""
    print("=" * 60)
    print("分子对接 - 受体制备")
    print("=" * 60)
    
    # 创建目录
    PDB_RAW_DIR.mkdir(parents=True, exist_ok=True)
    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\nPDB原始文件目录: {PDB_RAW_DIR}")
    print(f"受体输出目录: {RECEPTOR_DIR}")
    print(f"蛋白数: {len(RECEPTORS)}")
    
    # 1. 下载PDB文件
    print("\n[1/2] 下载PDB结构")
    pdb_files = {}
    for name, pdb_id in RECEPTORS.items():
        pdb_file = download_pdb(pdb_id, PDB_RAW_DIR)
        if pdb_file:
            pdb_files[name] = pdb_file
    
    print(f"\n成功下载 {len(pdb_files)}/{len(RECEPTORS)} 个PDB文件")
    
    # 2. 预处理受体
    print("\n[2/2] 预处理受体蛋白")
    receptor_files = {}
    for name, pdb_file in pdb_files.items():
        pdbqt_file = prepare_receptor(pdb_file, RECEPTOR_DIR, name)
        if pdbqt_file:
            receptor_files[name] = pdbqt_file
    
    # 总结
    print("\n" + "=" * 60)
    print("受体制备完成")
    print(f"  成功预处理: {len(receptor_files)}/{len(RECEPTORS)}")
    print(f"  输出目录: {RECEPTOR_DIR}")
    
    for name in RECEPTORS:
        status = "OK" if name in receptor_files else "FAILED"
        print(f"    {name:12s} [{status}]")
    
    print("\n下一步: 运行 02_prepare_ligand.py 准备配体文件")
    print("=" * 60)
    
    return 0 if len(receptor_files) > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
