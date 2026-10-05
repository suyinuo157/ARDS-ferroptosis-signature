# -*- coding: utf-8 -*-
"""
---
title: "分子对接：运行AutoDock Vina对接"
description: "使用AutoDock Vina对铁死亡关键蛋白与候选药物进行分子对接，计算结合亲和力（binding affinity），筛选高亲和力的药物-靶点对"
input: "results/15_docking/receptors/*.pdbqt, results/15_docking/ligands/*.pdbqt"
output: "results/15_docking/docking_results.csv, results/15_docking/poses/*.pdbqt"
dependencies: "vina (AutoDock Vina), pymol（可选，用于可视化）"
---

用法：
    python 03_run_docking.py
    python 03_run_docking.py --receptor GPX4 --ligand Sulforaphane

依赖：
    - AutoDock Vina (vina 或 vina_1.2.x)
    - OpenBabel (可选，用于分析)
"""

import os
import sys
import subprocess
import csv
import argparse
from pathlib import Path
from itertools import product


# ===================== 配置 =====================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

RECEPTOR_DIR = PROJECT_ROOT / "results" / "15_docking" / "receptors"
LIGAND_DIR = PROJECT_ROOT / "results" / "15_docking" / "ligands"
OUTPUT_DIR = PROJECT_ROOT / "results" / "15_docking" / "poses"
RESULT_CSV = PROJECT_ROOT / "results" / "15_docking" / "docking_results.csv"

# 每个蛋白的对接盒子参数（中心坐标 + 尺寸）
# 基于活性位点的大致位置，单位：Å
BOX_PARAMS = {
    "GPX4": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 25, "size_y": 25, "size_z": 25
    },
    "NFE2L2": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 30, "size_y": 30, "size_z": 30
    },
    "SLC7A11": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 25, "size_y": 25, "size_z": 25
    },
    "ACSL4": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 25, "size_y": 25, "size_z": 25
    },
    "HMOX1": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 25, "size_y": 25, "size_z": 25
    },
    "PTGS2": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 25, "size_y": 25, "size_z": 25
    },
    "TFRC": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 25, "size_y": 25, "size_z": 25
    },
    "HSPB1": {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 25, "size_y": 25, "size_z": 25
    }
}

# Vina可执行文件
VINA_CMD = os.environ.get("VINA_CMD", "vina")
VINA_CPU = 4
VINA_EXHAUSTIVENESS = 8
VINA_NUM_MODES = 9


def find_vina():
    """查找Vina可执行文件"""
    try:
        result = subprocess.run([VINA_CMD, "--version"], capture_output=True, text=True, timeout=10)
        print(f"  找到 Vina: {result.stdout.strip() or result.stderr.strip()}")
        return True
    except FileNotFoundError:
        print(f"  [错误] 未找到 Vina: {VINA_CMD}")
        print(f"  请安装 AutoDock Vina 并设置 VINA_CMD 环境变量")
        print(f"  下载地址: https://vina.scripps.edu/")
        return False
    except Exception as e:
        print(f"  [错误] Vina 调用失败: {e}")
        return False


def run_docking(receptor_file, ligand_file, output_file, box_params, log_file=None):
    """
    运行AutoDock Vina分子对接
    
    Args:
        receptor_file: 受体PDBQT文件路径
        ligand_file: 配体PDBQT文件路径
        output_file: 输出对接构象PDBQT文件路径
        box_params: 对接盒子参数字典
        log_file: 日志文件路径
    
    Returns:
        dict: 对接结果（结合能等）
    """
    # 构建命令
    cmd = [
        VINA_CMD,
        "--receptor", str(receptor_file),
        "--ligand", str(ligand_file),
        "--out", str(output_file),
        "--center_x", str(box_params["center_x"]),
        "--center_y", str(box_params["center_y"]),
        "--center_z", str(box_params["center_z"]),
        "--size_x", str(box_params["size_x"]),
        "--size_y", str(box_params["size_y"]),
        "--size_z", str(box_params["size_z"]),
        "--cpu", str(VINA_CPU),
        "--exhaustiveness", str(VINA_EXHAUSTIVENESS),
        "--num_modes", str(VINA_NUM_MODES)
    ]
    
    if log_file:
        cmd.extend(["--log", str(log_file)])
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        # 解析输出获取结合能
        binding_affinity = None
        for line in result.stdout.split("\n") + result.stderr.split("\n"):
            if "1" in line and "-" in line:
                parts = line.split()
                if len(parts) >= 4 and parts[0] == "1":
                    try:
                        binding_affinity = float(parts[1])
                    except ValueError:
                        pass
        
        return {
            "success": True,
            "binding_affinity": binding_affinity,
            "stdout": result.stdout,
            "stderr": result.stderr
        }
    
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "Timeout"}
    except Exception as e:
        return {"success": False, "error": str(e)}


def parse_vina_output(output_file):
    """解析Vina输出文件，获取所有构象的结合能"""
    if not output_file.exists():
        return []
    
    results = []
    with open(output_file) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        affinity = float(parts[3])
                        rmsd_lb = float(parts[4])
                        rmsd_ub = float(parts[5])
                        results.append({
                            "mode": len(results) + 1,
                            "affinity_kcal_mol": affinity,
                            "rmsd_lb": rmsd_lb,
                            "rmsd_ub": rmsd_ub
                        })
                    except (ValueError, IndexError):
                        pass
    
    return results


def main():
    """主函数：运行所有分子对接"""
    parser = argparse.ArgumentParser(description="AutoDock Vina 分子对接")
    parser.add_argument("--receptor", type=str, default=None,
                        help="指定单个受体蛋白名称")
    parser.add_argument("--ligand", type=str, default=None,
                        help="指定单个配体名称")
    args = parser.parse_args()
    
    print("=" * 60)
    print("分子对接 - AutoDock Vina")
    print("=" * 60)
    
    # 创建输出目录
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # 检查Vina
    print("\n检查Vina...")
    if not find_vina():
        print("\n[跳过] Vina未安装，无法执行对接")
        print("请安装 AutoDock Vina 后重新运行")
        print("=" * 60)
        return 1
    
    # 收集受体和配体
    receptors = {}
    if args.receptor:
        rec_file = RECEPTOR_DIR / f"{args.receptor}.pdbqt"
        if rec_file.exists():
            receptors[args.receptor] = rec_file
    else:
        for f in RECEPTOR_DIR.glob("*.pdbqt"):
            receptors[f.stem] = f
    
    ligands = {}
    if args.ligand:
        lig_file = LIGAND_DIR / f"{args.ligand}.pdbqt"
        if lig_file.exists():
            ligands[args.ligand] = lig_file
    else:
        for f in LIGAND_DIR.glob("*.pdbqt"):
            ligands[f.stem] = f
    
    print(f"\n受体数: {len(receptors)}")
    print(f"配体数: {len(ligands)}")
    print(f"总对接数: {len(receptors) * len(ligands)}")
    
    if len(receptors) == 0 or len(ligands) == 0:
        print("\n[错误] 受体或配体文件为空")
        print(f"  受体目录: {RECEPTOR_DIR}")
        print(f"  配体目录: {LIGAND_DIR}")
        return 1
    
    # 运行对接
    print("\n" + "-" * 60)
    print("运行分子对接...")
    print("-" * 60)
    
    all_results = []
    total = len(receptors) * len(ligands)
    done = 0
    
    for rec_name, rec_file in receptors.items():
        box = BOX_PARAMS.get(rec_name, BOX_PARAMS["GPX4"])
        
        for lig_name, lig_file in ligands.items():
            done += 1
            print(f"\n[{done}/{total}] {rec_name} + {lig_name}")
            
            # 输出文件
            out_name = f"{rec_name}_{lig_name}"
            out_file = OUTPUT_DIR / f"{out_name}.pdbqt"
            log_file = OUTPUT_DIR / f"{out_name}.log"
            
            # 跳过已完成的
            if out_file.exists() and out_file.stat().st_size > 100:
                # 解析已有结果
                poses = parse_vina_output(out_file)
                if poses:
                    best = poses[0]
                    all_results.append({
                        "receptor": rec_name,
                        "ligand": lig_name,
                        "best_affinity_kcal_mol": best["affinity_kcal_mol"],
                        "best_rmsd_lb": best["rmsd_lb"],
                        "best_rmsd_ub": best["rmsd_ub"],
                        "num_modes": len(poses),
                        "status": "done"
                    })
                    print(f"  [skip] 已存在, best affinity = {best['affinity_kcal_mol']} kcal/mol")
                    continue
            
            # 运行对接
            result = run_docking(rec_file, lig_file, out_file, box, log_file)
            
            if result["success"]:
                poses = parse_vina_output(out_file)
                if poses:
                    best = poses[0]
                    all_results.append({
                        "receptor": rec_name,
                        "ligand": lig_name,
                        "best_affinity_kcal_mol": best["affinity_kcal_mol"],
                        "best_rmsd_lb": best["rmsd_lb"],
                        "best_rmsd_ub": best["rmsd_ub"],
                        "num_modes": len(poses),
                        "status": "success"
                    })
                    print(f"  完成, best affinity = {best['affinity_kcal_mol']} kcal/mol")
                else:
                    all_results.append({
                        "receptor": rec_name,
                        "ligand": lig_name,
                        "best_affinity_kcal_mol": "",
                        "best_rmsd_lb": "",
                        "best_rmsd_ub": "",
                        "num_modes": 0,
                        "status": "no_result"
                    })
                    print(f"  完成但未解析到结果")
            else:
                all_results.append({
                    "receptor": rec_name,
                    "ligand": lig_name,
                    "best_affinity_kcal_mol": "",
                    "best_rmsd_lb": "",
                    "best_rmsd_ub": "",
                    "num_modes": 0,
                    "status": f"failed: {result.get('error', 'unknown')}"
                })
                print(f"  失败: {result.get('error', 'unknown')}")
    
    # 保存结果
    print("\n" + "=" * 60)
    print("保存结果...")
    
    if all_results:
        with open(RESULT_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
            writer.writeheader()
            writer.writerows(all_results)
        
        print(f"  结果已保存: {RESULT_CSV}")
    
    # 统计
    success_count = sum(1 for r in all_results if r["status"] == "success" or r["status"] == "done")
    print(f"\n对接完成: {success_count}/{total} 成功")
    
    # Top 结果
    successful = [r for r in all_results if r["best_affinity_kcal_mol"]]
    if successful:
        successful.sort(key=lambda x: float(x["best_affinity_kcal_mol"]))
        print("\nTop 10 结合亲和力:")
        for i, r in enumerate(successful[:10]):
            print(f"  {i+1:2d}. {r['receptor']:10s} + {r['ligand']:25s} = "
                  f"{r['best_affinity_kcal_mol']:>6s} kcal/mol")
    
    print("\n" + "=" * 60)
    print("分子对接完成")
    print(f"  结果文件: {RESULT_CSV}")
    print(f"  构象目录: {OUTPUT_DIR}")
    print("=" * 60)
    
    return 0 if success_count > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
