# -*- coding: utf-8 -*-
"""
---
title: "分子对接：配体制备（小分子配体预处理）"
description: "从PubChem下载候选药物的3D结构（SDF格式），使用OpenBabel或AutoDockTools转换为PDBQT格式，添加Gasteiger电荷和可旋转键，准备用于分子对接的配体文件"
input: "候选药物名称列表（Sulforaphane, Dimethyl fumarate, Ferrostatin-1, Erastin, Deferoxamine等）"
output: "results/15_docking/ligands/*.pdbqt （预处理后的配体文件）"
dependencies: "requests, openbabel (可选，用于格式转换)"
---

用法：
    python 02_prepare_ligand.py

依赖：
    - OpenBabel (用于格式转换) 或 MGLTools prepare_ligand4.py
    - requests (用于PubChem下载)
"""

import os
import sys
import subprocess
import requests
from pathlib import Path


# ===================== 配置 =====================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

LIGAND_DIR = PROJECT_ROOT / "results" / "15_docking" / "ligands"
SDF_RAW_DIR = PROJECT_ROOT / "results" / "15_docking" / "sdf_raw"

# 候选药物 - PubChem CID 映射
LIGANDS = {
    # Nrf2 激活剂
    "Sulforaphane": "5350",
    "Dimethyl_fumarate": "6375",
    "Oltipraz": "3757",
    "Resveratrol": "445154",
    "CDDO_Methyl_ester": "10032964",
    
    # 铁死亡抑制剂（GPX4相关）
    "Ferrostatin-1": "5281936",
    "Liproxstatin-1": "118702733",
    "RSL3": "11725265",
    
    # SLC7A11 抑制剂
    "Erastin": "65399",
    "Sulfasalazine": "5395",
    "Sorafenib": "216239",
    
    # 铁螯合剂
    "Deferoxamine": "4261",
    "Deferasirox": "3731",
    
    # ACSL4 抑制剂
    "Rosiglitazone": "77999",
    "Triacsin_C": "93442",
    
    # 抗炎药
    "Dexamethasone": "5743",
    "Methylprednisolone": "6741",
    
    # 抗氧化剂
    "N-acetylcysteine": "12015",
    "Vitamin_E": "14985",
    "Vitamin_C": "54670067",
    
    # HMOX1 诱导剂
    "Hemin": "5362",
    "Curcumin": "969516"
}


def download_sdf(drug_name, cid, output_dir):
    """从PubChem下载SDF文件"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
    output_file = output_dir / f"{drug_name}.sdf"
    
    if output_file.exists() and output_file.stat().st_size > 100:
        print(f"  [skip] {drug_name}.sdf 已存在")
        return output_file
    
    print(f"  下载 {drug_name} (CID:{cid}) ...", end=" ")
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        output_file.write_bytes(response.content)
        print("完成")
        return output_file
    except Exception as e:
        print(f"失败: {e}")
        return None


def sdf_to_pdbqt(sdf_file, output_dir, ligand_name):
    """
    使用OpenBabel将SDF转换为PDBQT格式
    - 添加Gasteiger电荷
    - 识别可旋转键
    """
    pdbqt_file = output_dir / f"{ligand_name}.pdbqt"
    
    if pdbqt_file.exists():
        print(f"  [skip] {ligand_name}.pdbqt 已存在")
        return pdbqt_file
    
    # 尝试使用openbabel
    try:
        cmd = [
            "obabel",
            "-isdf", str(sdf_file),
            "-opdbqt",
            "-O", str(pdbqt_file),
            "-h",           # 添加氢
            "-p", "7.4",    # pH 7.4 质子化
            "--partialcharge", "gasteiger"  # Gasteiger电荷
        ]
        
        print(f"  OpenBabel转换 {ligand_name} ...", end=" ")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if pdbqt_file.exists() and pdbqt_file.stat().st_size > 100:
            print("完成")
            return pdbqt_file
        else:
            print("失败")
            if result.stderr:
                print(f"  stderr: {result.stderr.strip()[:100]}")
    except FileNotFoundError:
        print("  [warn] OpenBabel未安装")
    except Exception as e:
        print(f"  错误: {e}")
    
    # 备用：尝试使用MGLTools
    mgl_path = os.environ.get("MGLTOOLS_PATH", "/usr/local/MGLToolsPckgs")
    prep_ligand = os.path.join(mgl_path, "AutoDockTools/Utilities24/prepare_ligand4.py")
    
    if os.path.exists(prep_ligand):
        # 先转PDB
        pdb_file = output_dir / f"{ligand_name}.pdb"
        try:
            cmd_pdb = ["obabel", "-isdf", str(sdf_file), "-opdb", "-O", str(pdb_file), "-h"]
            subprocess.run(cmd_pdb, capture_output=True, timeout=60)
        except:
            pass
        
        if pdb_file.exists():
            cmd = [
                "pythonsh", prep_ligand,
                "-l", str(pdb_file),
                "-o", str(pdbqt_file),
                "-A", "hydrogens",
                "-v"
            ]
            try:
                print(f"  MGLTools转换 {ligand_name} ...", end=" ")
                subprocess.run(cmd, capture_output=True, timeout=120)
                if pdbqt_file.exists():
                    print("完成")
                    return pdbqt_file
                else:
                    print("失败")
            except Exception as e:
                print(f"错误: {e}")
    
    print(f"  [note] 请手动使用AutoDockTools或OpenBabel处理 {sdf_file}")
    print(f"         输出至: {pdbqt_file}")
    return None


def main():
    """主函数：下载并预处理所有配体"""
    print("=" * 60)
    print("分子对接 - 配体制备")
    print("=" * 60)
    
    SDF_RAW_DIR.mkdir(parents=True, exist_ok=True)
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\nSDF原始文件目录: {SDF_RAW_DIR}")
    print(f"配体输出目录: {LIGAND_DIR}")
    print(f"配体数: {len(LIGANDS)}")
    
    # 1. 下载SDF文件
    print("\n[1/2] 从PubChem下载SDF结构")
    sdf_files = {}
    for name, cid in LIGANDS.items():
        sdf_file = download_sdf(name, cid, SDF_RAW_DIR)
        if sdf_file:
            sdf_files[name] = sdf_file
    
    print(f"\n成功下载 {len(sdf_files)}/{len(LIGANDS)} 个SDF文件")
    
    # 2. 转换为PDBQT
    print("\n[2/2] SDF -> PDBQT 转换")
    ligand_files = {}
    for name, sdf_file in sdf_files.items():
        pdbqt_file = sdf_to_pdbqt(sdf_file, LIGAND_DIR, name)
        if pdbqt_file:
            ligand_files[name] = pdbqt_file
    
    # 总结
    print("\n" + "=" * 60)
    print("配体制备完成")
    print(f"  成功预处理: {len(ligand_files)}/{len(LIGANDS)}")
    print(f"  输出目录: {LIGAND_DIR}")
    
    print("\n配体列表:")
    for name in LIGANDS:
        status = "OK" if name in ligand_files else "FAILED"
        print(f"    {name:30s} [{status}]")
    
    print("\n下一步: 运行 03_run_docking.py 执行分子对接")
    print("=" * 60)
    
    return 0 if len(ligand_files) > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
