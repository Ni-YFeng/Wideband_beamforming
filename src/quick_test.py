"""
快速测试脚本 - 使用简化参数
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys, os

# 设置中文字体
rcParams['font.sans-serif'] = ['SimHei']
rcParams['axes.unicode_minus'] = False

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from system_model import CellFreeMIMOSystem
from precoding_algorithm import DistributedPrecoder, BaselineSchemes

def quick_test():
    """快速测试基本功能"""
    print("=" * 60)
    print("快速测试 - Cell-Free MIMO分布式预编码系统")
    print("=" * 60)

    # 简化配置
    config = {
        'B': 2,  # 减少到2个AP
        'K': 2,  # 减少到2个用户
        'Nt': 64,  # 减少天线数
        'NRF': 2,  # 减少RF链
        'KD': 4,   # 减少TTD组
        'M': 16,   # 大幅减少子载波
        'fc': 300e9,
        'W': 30e9
    }

    # 创建系统
    system = CellFreeMIMOSystem(config)

    # 生成信道
    print("\n生成信道...")
    H, _ = system.generate_channel_cluster_based(Lc=2, Lp=5)
    print(f"信道形状: {H.shape}")

    # 测试预编码器初始化
    print("\n测试预编码器初始化...")
    precoder = DistributedPrecoder(system)
    A, T, D = precoder.initialize_precoder(H)
    print(f"相移矩阵A[0]形状: {A[0].shape}")
    print(f"TTD矩阵T[0]形状: {T[0].shape}")
    print(f"数字预编码D[0]形状: {D[0].shape}")

    # 运行算法（仅5次迭代）
    print("\n运行DP-AltMin算法 (5次迭代)...")
    results = precoder.distributed_dp_altmin(H, max_iter=5)

    print(f"\n收敛和速率序列: {[f'{sr:.2f}' for sr in results['sum_rates']]}")

    # 计算性能
    SINR = system.compute_sinr(H, results['A'], results['T'], results['D'])
    sum_rate = system.compute_sum_rate(SINR)

    print(f"\n最终和速率: {sum_rate:.2f} bps/Hz")

    # 测试基准方案
    print("\n测试基准方案...")
    baselines = BaselineSchemes(system)

    print("\n全数字预编码...")
    res_fd = baselines.fully_digital(H)
    print(f"全数字和速率: {res_fd['sum_rate']:.2f} bps/Hz")

    print("\n传统混合预编码...")
    res_th = baselines.traditional_hybrid(H)
    print(f"传统混合和速率: {res_th['sum_rate']:.2f} bps/Hz")

    print("\n非协调本地预编码...")
    res_uc = baselines.uncoordinated_local(H)
    print(f"非协调和速率: {res_uc['sum_rate']:.2f} bps/Hz")

    print("\n" + "=" * 60)
    print("测试完成！")
    print("=" * 60)

    # 绘制简单的收敛图
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(results['sum_rates'])+1), results['sum_rates'],
             'b-o', linewidth=2, markersize=6)
    plt.xlabel('迭代次数', fontsize=12)
    plt.ylabel('和速率 (bps/Hz)', fontsize=12)
    plt.title('DP-AltMin算法收敛曲线 (快速测试)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('quick_test_convergence.png', dpi=300, bbox_inches='tight')
    print("\n收敛曲线已保存为: quick_test_convergence.png")

    # 性能对比柱状图
    schemes = ['DP-AltMin', '全数字', '传统混合', '非协调']
    rates = [sum_rate, res_fd['sum_rate'], res_th['sum_rate'], res_uc['sum_rate']]

    plt.figure(figsize=(8, 5))
    bars = plt.bar(schemes, rates, color=['blue', 'red', 'green', 'purple'],
                   alpha=0.7, edgecolor='black')
    plt.ylabel('和速率 (bps/Hz)', fontsize=12)
    plt.title('不同方案的性能对比 (快速测试)', fontsize=14)
    plt.grid(True, alpha=0.3, axis='y')

    # 在柱子上添加数值标签
    for bar, rate in zip(bars, rates):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{rate:.2f}', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig('quick_test_comparison.png', dpi=300, bbox_inches='tight')
    print("性能对比图已保存为: quick_test_comparison.png")

    return results

if __name__ == "__main__":
    results = quick_test()
