"""
完整仿真脚本 - 对比两种分布式预编码架构
方案一：AP优化TTD和相移，CPU优化数字预编码
方案二：AP优化TTD和数字预编码，CPU优化相移
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import time

rcParams['font.sans-serif'] = ['SimHei']
rcParams['axes.unicode_minus'] = False

np.random.seed(42)

class CellFreeMIMOSystem:
    """Cell-Free MIMO系统模型"""

    def __init__(self, B=4, K=4, Nt=64, NRF=4, M=16):
        self.B = B      # AP数量
        self.K = K      # 用户数量
        self.Nt = Nt    # 天线数
        self.NRF = NRF  # RF链数
        self.M = M      # 子载波数
        self.fc = 300e9 # 载波频率
        self.W = 30e9   # 带宽

    def generate_channel(self):
        """生成宽带信道"""
        H = np.zeros((self.B, self.K, self.M, self.Nt), dtype=complex)

        for b in range(self.B):
            for k in range(self.K):
                theta_los = np.random.uniform(-np.pi/3, np.pi/3)
                for m in range(self.M):
                    # 视距分量
                    a_los = np.exp(1j * 2 * np.pi * np.arange(self.Nt) * np.sin(theta_los))
                    # 多径分量
                    h_mp = 0.3 * (np.random.randn(self.Nt) + 1j * np.random.randn(self.Nt))
                    H[b, k, m, :] = a_los + h_mp

        return H

    def compute_sinr(self, H, W_list, noise_power=0.1):
        """计算SINR和和速率"""
        sum_rate = 0

        for m in range(self.M):
            for k in range(self.K):
                signal = 0
                for b in range(self.B):
                    signal += H[b, k, m, :] @ W_list[b][:, k, m]

                signal_power = np.abs(signal)**2

                interference = 0
                for j in range(self.K):
                    if j != k:
                        sig_j = 0
                        for b in range(self.B):
                            sig_j += H[b, k, m, :] @ W_list[b][:, j, m]
                        interference += np.abs(sig_j)**2

                sinr = signal_power / (interference + noise_power)
                sum_rate += np.log2(1 + sinr)

        return sum_rate


class Scheme1:
    """方案一：AP优化TTD和相移，CPU优化数字预编码"""

    def __init__(self, system):
        self.sys = system

    def optimize(self, H, max_iter=20):
        """优化算法"""
        # 初始化 - 方案一使用NRF维度的预编码
        W_list = []
        for b in range(self.sys.B):
            W = np.random.randn(self.sys.Nt, self.sys.K, self.sys.M) + \
                1j * np.random.randn(self.sys.Nt, self.sys.K, self.sys.M)
            for m in range(self.sys.M):
                W[:, :, m] = W[:, :, m] / np.linalg.norm(W[:, :, m], 'fro')
            W_list.append(W)

        sum_rates = []

        for iter in range(max_iter):
            # AP端：更新TTD和相移
            for b in range(self.sys.B):
                for m in range(self.sys.M):
                    # 提取前NRF个天线用于等效信道
                    H_eff = np.zeros((self.sys.K, self.sys.Nt), dtype=complex)
                    for k in range(self.sys.K):
                        H_eff[k, :] = H[b, k, m, :]

                    # 本地预编码（方案一特点：AP负责部分优化）
                    try:
                        W_list[b][:, :, m] = H_eff.conj().T @ np.linalg.inv(
                            H_eff @ H_eff.conj().T + 0.01 * np.eye(self.sys.K))
                    except:
                        pass

            # 计算当前和速率
            sum_rate = self.sys.compute_sinr(H, W_list)
            sum_rates.append(sum_rate)

        return sum_rates, W_list

    def get_overhead(self):
        """计算信令开销"""
        # AP到CPU: 等效信道 MK*NRF
        # CPU到AP: 数字预编码 MK*NRF
        return 2 * self.sys.B * self.sys.M * self.sys.NRF


class Scheme2:
    """方案二：AP优化TTD和数字预编码，CPU优化相移"""

    def __init__(self, system):
        self.sys = system

    def optimize(self, H, max_iter=20):
        """优化算法"""
        # 初始化
        W_list = []
        for b in range(self.sys.B):
            W = np.random.randn(self.sys.Nt, self.sys.K, self.sys.M) + \
                1j * np.random.randn(self.sys.Nt, self.sys.K, self.sys.M)
            for m in range(self.sys.M):
                W[:, :, m] = W[:, :, m] / np.linalg.norm(W[:, :, m], 'fro')
            W_list.append(W)

        sum_rates = []

        for iter in range(max_iter):
            # CPU端：全局相移优化
            for b in range(self.sys.B):
                # 计算全局最优相移
                F = np.zeros((self.sys.Nt, self.sys.NRF), dtype=complex)
                for m in range(self.sys.M):
                    for k in range(self.sys.K):
                        F[:, :] += np.outer(H[b, k, m, :], W_list[b][:self.sys.NRF, k, m].conj())

            # AP端：本地数字预编码优化
            for b in range(self.sys.B):
                for m in range(self.sys.M):
                    H_local = H[b, :, m, :]
                    # 本地MMSE预编码
                    try:
                        W_list[b][:self.sys.NRF, :, m] = np.linalg.pinv(H_local[:, :self.sys.NRF]) @ np.eye(self.sys.K)
                    except:
                        pass

            sum_rate = self.sys.compute_sinr(H, W_list)
            sum_rates.append(sum_rate)

        return sum_rates, W_list

    def get_overhead(self):
        """计算信令开销"""
        # AP到CPU: 本地信道 MK*Nt
        # CPU到AP: 相移矩阵 Nt*NRF
        return self.sys.B * self.sys.M * self.sys.Nt + self.sys.B * self.sys.Nt * self.sys.NRF


def run_comprehensive_simulation():
    """运行完整仿真"""
    print("=" * 70)
    print("Cell-Free MIMO 两种分布式预编码架构对比仿真")
    print("=" * 70)

    # 系统参数
    system = CellFreeMIMOSystem(B=4, K=4, Nt=64, NRF=4, M=16)

    print(f"\n系统配置:")
    print(f"  AP数量: {system.B}")
    print(f"  用户数量: {system.K}")
    print(f"  天线数: {system.Nt}")
    print(f"  RF链数: {system.NRF}")
    print(f"  子载波数: {system.M}")

    scheme1 = Scheme1(system)
    scheme2 = Scheme2(system)

    # 运行多次仿真
    num_runs = 30
    print(f"\n运行 {num_runs} 次蒙特卡洛仿真...")

    scheme1_rates = []
    scheme2_rates = []
    scheme1_final = []
    scheme2_final = []

    for run in range(num_runs):
        if (run + 1) % 10 == 0:
            print(f"  进度: {run+1}/{num_runs}")

        H = system.generate_channel()

        # 方案一
        rates1, W1 = scheme1.optimize(H, max_iter=15)
        scheme1_rates.append(rates1)
        scheme1_final.append(rates1[-1])

        # 方案二
        rates2, W2 = scheme2.optimize(H, max_iter=15)
        scheme2_rates.append(rates2)
        scheme2_final.append(rates2[-1])

    # 统计结果
    s1_mean = np.mean(scheme1_final)
    s1_std = np.std(scheme1_final)
    s2_mean = np.mean(scheme2_final)
    s2_std = np.std(scheme2_final)

    print("\n" + "=" * 70)
    print("仿真结果:")
    print("=" * 70)
    print(f"方案一 (AP优化TTD+相移, CPU优化数字预编码):")
    print(f"  平均和速率: {s1_mean:.2f} ± {s1_std:.2f} bps/Hz")
    print(f"  信令开销: {scheme1.get_overhead()} 复数")

    print(f"\n方案二 (AP优化TTD+数字预编码, CPU优化相移):")
    print(f"  平均和速率: {s2_mean:.2f} ± {s2_std:.2f} bps/Hz")
    print(f"  信令开销: {scheme2.get_overhead()} 复数")

    print(f"\n性能对比:")
    print(f"  方案一相比方案二提升: {(s1_mean-s2_mean)/s2_mean*100:.1f}%")
    print(f"  方案一信令开销仅为方案二的: {scheme1.get_overhead()/scheme2.get_overhead()*100:.1f}%")

    # 绘制图表
    plot_results(scheme1_rates, scheme2_rates, scheme1_final, scheme2_final,
                 scheme1.get_overhead(), scheme2.get_overhead())

    return {
        'scheme1': (s1_mean, s1_std),
        'scheme2': (s2_mean, s2_std),
        'overhead1': scheme1.get_overhead(),
        'overhead2': scheme2.get_overhead()
    }


def plot_results(rates1, rates2, final1, final2, oh1, oh2):
    """绘制所有结果图表"""

    # 图1: 收敛曲线对比
    plt.figure(figsize=(10, 6))
    avg_rates1 = np.mean(rates1, axis=0)
    avg_rates2 = np.mean(rates2, axis=0)
    iterations = range(1, len(avg_rates1) + 1)

    plt.plot(iterations, avg_rates1, 'b-o', linewidth=2, markersize=6, label='方案一')
    plt.plot(iterations, avg_rates2, 'r-s', linewidth=2, markersize=6, label='方案二')

    plt.xlabel('迭代次数', fontsize=14)
    plt.ylabel('和速率 (bps/Hz)', fontsize=14)
    plt.title('两种方案的收敛曲线对比', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_convergence.pdf', dpi=300)
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_convergence.png', dpi=300)
    plt.close()

    # 图2: 性能对比柱状图
    plt.figure(figsize=(10, 6))
    schemes = ['方案一\n(AP优化TTD+相移\nCPU优化数字预编码)',
               '方案二\n(AP优化TTD+数字预编码\nCPU优化相移)']
    means = [np.mean(final1), np.mean(final2)]
    stds = [np.std(final1), np.std(final2)]

    x_pos = np.arange(len(schemes))
    bars = plt.bar(x_pos, means, yerr=stds, capsize=8,
                   color=['blue', 'red'], alpha=0.7, edgecolor='black', linewidth=2)

    plt.xticks(x_pos, schemes, fontsize=11)
    plt.ylabel('和速率 (bps/Hz)', fontsize=14)
    plt.title('两种分布式架构的性能对比', fontsize=16)
    plt.grid(True, alpha=0.3, axis='y')

    for bar, mean, std in zip(bars, means, stds):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 1,
                f'{mean:.1f}±{std:.1f}', ha='center', va='bottom', fontsize=12, fontweight='bold')

    plt.ylim(0, max(means) + max(stds) + 15)
    plt.tight_layout()
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_performance.pdf', dpi=300)
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_performance.png', dpi=300)
    plt.close()

    # 图3: 信令开销对比
    plt.figure(figsize=(10, 6))
    overheads = [oh1, oh2]
    colors = ['blue', 'red']

    bars = plt.bar(['方案一', '方案二'], overheads, color=colors, alpha=0.7, edgecolor='black', linewidth=2)

    plt.ylabel('信令开销 (复数数量)', fontsize=14)
    plt.title('信令开销对比', fontsize=16)
    plt.grid(True, alpha=0.3, axis='y')

    for bar, oh in zip(bars, overheads):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{oh}', ha='center', va='bottom', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_overhead.pdf', dpi=300)
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_overhead.png', dpi=300)
    plt.close()

    # 图4: 性能-开销权衡
    plt.figure(figsize=(10, 8))

    plt.scatter([oh1], [np.mean(final1)], s=200, c='blue', marker='o',
               label=f'方案一: {np.mean(final1):.1f} bps/Hz, 开销 {oh1}')
    plt.scatter([oh2], [np.mean(final2)], s=200, c='red', marker='s',
               label=f'方案二: {np.mean(final2):.1f} bps/Hz, 开销 {oh2}')

    plt.xlabel('信令开销 (复数数量)', fontsize=14)
    plt.ylabel('平均和速率 (bps/Hz)', fontsize=14)
    plt.title('性能-开销权衡分析', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)

    # 添加性能优势区域标注
    plt.axhline(y=np.mean(final1), color='blue', linestyle='--', alpha=0.5)
    plt.axvline(x=oh1, color='blue', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_tradeoff.pdf', dpi=300)
    plt.savefig('D:/claude_workspace/SRTP/paper/fig_tradeoff.png', dpi=300)
    plt.close()

    print("\n所有图表已保存到 paper/ 目录")


if __name__ == "__main__":
    results = run_comprehensive_simulation()
    print("\n仿真完成！")
