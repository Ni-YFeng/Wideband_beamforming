"""
完整的仿真实验脚本
生成论文所需的所有性能对比图表
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import json
import time
import sys
import os

# 设置中文字体
rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 导入自定义模块
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from system_model import CellFreeMIMOSystem
from precoding_algorithm import DistributedPrecoder, BaselineSchemes

class SimulationExperiment:
    """仿真实验类"""

    def __init__(self, output_dir: str = "../results"):
        """
        初始化仿真实验

        Parameters:
        -----------
        output_dir : str
            结果输出目录
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        # 默认系统配置
        self.config = {
            'B': 4,
            'K': 4,
            'Nt': 256,
            'NRF': 4,
            'KD': 16,
            'M': 128,
            'fc': 300e9,
            'W': 30e9,
            'Pmax': 1.0
        }

        print("=" * 80)
        print("Cell-Free MIMO 分布式预编码与时延-相位联合优化仿真系统")
        print("=" * 80)
        print(f"系统配置:")
        for key, value in self.config.items():
            print(f"  {key}: {value}")
        print("=" * 80)

    def run_convergence_analysis(self, num_runs: int = 5):
        """
        运行收敛性分析

        Parameters:
        -----------
        num_runs : int
            独立运行次数
        """
        print("\n" + "=" * 80)
        print("实验1: 算法收敛性分析")
        print("=" * 80)

        system = CellFreeMIMOSystem(self.config)
        precoder = DistributedPrecoder(system)

        all_sum_rates = []

        for run in range(num_runs):
            print(f"\n运行第 {run+1}/{num_runs} 次...")
            H, _ = system.generate_channel_cluster_based()
            results = precoder.distributed_dp_altmin(H, max_iter=50)
            all_sum_rates.append(results['sum_rates'])

        # 计算平均收敛曲线
        max_len = max(len(sr) for sr in all_sum_rates)
        avg_sum_rates = np.zeros(max_len)
        for i in range(max_len):
            values = []
            for sr in all_sum_rates:
                if i < len(sr):
                    values.append(sr[i])
                else:
                    values.append(sr[-1])
            avg_sum_rates[i] = np.mean(values)

        # 绘制收敛曲线
        plt.figure(figsize=(10, 6))
        iterations = range(1, len(avg_sum_rates) + 1)
        plt.plot(iterations, avg_sum_rates, 'b-o', linewidth=2, markersize=4)

        # 添加标准差区域
        std_sum_rates = []
        for i in range(max_len):
            values = []
            for sr in all_sum_rates:
                if i < len(sr):
                    values.append(sr[i])
                else:
                    values.append(sr[-1])
            std_sum_rates.append(np.std(values))

        plt.fill_between(iterations,
                        avg_sum_rates - std_sum_rates,
                        avg_sum_rates + std_sum_rates,
                        alpha=0.3, color='blue')

        plt.xlabel('迭代次数', fontsize=14)
        plt.ylabel('和速率 (bps/Hz)', fontsize=14)
        plt.title('DP-AltMin算法收敛曲线', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'convergence_curve.pdf'),
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir, 'convergence_curve.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"\n收敛曲线已保存")
        print(f"最终平均和速率: {avg_sum_rates[-1]:.2f} ± {std_sum_rates[-1]:.2f} bps/Hz")

        # 保存数据
        np.savez(os.path.join(self.output_dir, 'convergence_data.npz'),
                iterations=iterations,
                avg_sum_rates=avg_sum_rates,
                std_sum_rates=std_sum_rates)

    def run_snr_sweep(self, snr_range: np.ndarray = None, num_runs: int = 10):
        """
        运行SNR扫描实验

        Parameters:
        -----------
        snr_range : np.ndarray
            SNR范围 (dB)
        num_runs : int
            每个SNR点的独立运行次数
        """
        print("\n" + "=" * 80)
        print("实验2: 和速率 vs SNR")
        print("=" * 80)

        if snr_range is None:
            snr_range = np.arange(-10, 21, 5)  # -10 to 20 dB

        system = CellFreeMIMOSystem(self.config)
        precoder = DistributedPrecoder(system)
        baselines = BaselineSchemes(system)

        results_dp_altmin = []
        results_fully_digital = []
        results_traditional = []
        results_uncoordinated = []

        for snr_db in snr_range:
            print(f"\nSNR = {snr_db} dB")

            # 设置噪声功率
            snr_linear = 10 ** (snr_db / 10)
            noise_power = self.config['Pmax'] / snr_linear

            sum_rate_dp = []
            sum_rate_fd = []
            sum_rate_th = []
            sum_rate_uc = []

            for run in range(num_runs):
                H, _ = system.generate_channel_cluster_based()

                # DP-AltMin
                res_dp = precoder.distributed_dp_altmin(H, max_iter=30)
                sum_rate_dp.append(res_dp['final_sum_rate'])

                # 全数字
                res_fd = baselines.fully_digital(H)
                sum_rate_fd.append(res_fd['sum_rate'])

                # 传统混合
                res_th = baselines.traditional_hybrid(H)
                sum_rate_th.append(res_th['sum_rate'])

                # 非协调
                res_uc = baselines.uncoordinated_local(H)
                sum_rate_uc.append(res_uc['sum_rate'])

            results_dp_altmin.append(np.mean(sum_rate_dp))
            results_fully_digital.append(np.mean(sum_rate_fd))
            results_traditional.append(np.mean(sum_rate_th))
            results_uncoordinated.append(np.mean(sum_rate_uc))

            print(f"  DP-AltMin: {results_dp_altmin[-1]:.2f} bps/Hz")
            print(f"  全数字:   {results_fully_digital[-1]:.2f} bps/Hz")
            print(f"  传统混合: {results_traditional[-1]:.2f} bps/Hz")
            print(f"  非协调:   {results_uncoordinated[-1]:.2f} bps/Hz")

        # 绘制性能对比图
        plt.figure(figsize=(10, 6))
        plt.plot(snr_range, results_dp_altmin, 'b-s', linewidth=2,
                markersize=8, label='DP-AltMin (本文方案)')
        plt.plot(snr_range, results_fully_digital, 'r-^', linewidth=2,
                markersize=8, label='全数字预编码')
        plt.plot(snr_range, results_traditional, 'g-d', linewidth=2,
                markersize=8, label='传统混合预编码')
        plt.plot(snr_range, results_uncoordinated, 'm-o', linewidth=2,
                markersize=8, label='非协调本地预编码')

        plt.xlabel('SNR (dB)', fontsize=14)
        plt.ylabel('和速率 (bps/Hz)', fontsize=14)
        plt.title('不同SNR下的和速率性能对比', fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'sum_rate_vs_snr.pdf'),
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir, 'sum_rate_vs_snr.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"\n性能对比图已保存")

        # 保存数据
        np.savez(os.path.join(self.output_dir, 'snr_sweep_data.npz'),
                snr_range=snr_range,
                dp_altmin=results_dp_altmin,
                fully_digital=results_fully_digital,
                traditional=results_traditional,
                uncoordinated=results_uncoordinated)

    def run_user_sweep(self, user_range: list = None, snr_db: float = 0,
                       num_runs: int = 10):
        """
        运行用户数扫描实验

        Parameters:
        -----------
        user_range : list
            用户数范围
        snr_db : float
            固定SNR (dB)
        num_runs : int
            独立运行次数
        """
        print("\n" + "=" * 80)
        print("实验3: 和速率 vs 用户数")
        print("=" * 80)

        if user_range is None:
            user_range = [2, 4, 6, 8, 10]

        snr_linear = 10 ** (snr_db / 10)
        noise_power = self.config['Pmax'] / snr_linear

        results_dp_altmin = []
        results_fully_digital = []
        results_traditional = []
        results_uncoordinated = []

        for K in user_range:
            print(f"\n用户数 K = {K}")

            # 更新配置
            config = self.config.copy()
            config['K'] = K
            system = CellFreeMIMOSystem(config)
            precoder = DistributedPrecoder(system)
            baselines = BaselineSchemes(system)

            sum_rate_dp = []
            sum_rate_fd = []
            sum_rate_th = []
            sum_rate_uc = []

            for run in range(num_runs):
                H, _ = system.generate_channel_cluster_based()

                # DP-AltMin
                res_dp = precoder.distributed_dp_altmin(H, max_iter=30)
                sum_rate_dp.append(res_dp['final_sum_rate'])

                # 全数字
                res_fd = baselines.fully_digital(H)
                sum_rate_fd.append(res_fd['sum_rate'])

                # 传统混合
                res_th = baselines.traditional_hybrid(H)
                sum_rate_th.append(res_th['sum_rate'])

                # 非协调
                res_uc = baselines.uncoordinated_local(H)
                sum_rate_uc.append(res_uc['sum_rate'])

            results_dp_altmin.append(np.mean(sum_rate_dp))
            results_fully_digital.append(np.mean(sum_rate_fd))
            results_traditional.append(np.mean(sum_rate_th))
            results_uncoordinated.append(np.mean(sum_rate_uc))

            print(f"  DP-AltMin: {results_dp_altmin[-1]:.2f} bps/Hz")
            print(f"  全数字:   {results_fully_digital[-1]:.2f} bps/Hz")
            print(f"  传统混合: {results_traditional[-1]:.2f} bps/Hz")
            print(f"  非协调:   {results_uncoordinated[-1]:.2f} bps/Hz")

        # 绘制性能对比图
        plt.figure(figsize=(10, 6))
        plt.plot(user_range, results_dp_altmin, 'b-s', linewidth=2,
                markersize=8, label='DP-AltMin (本文方案)')
        plt.plot(user_range, results_fully_digital, 'r-^', linewidth=2,
                markersize=8, label='全数字预编码')
        plt.plot(user_range, results_traditional, 'g-d', linewidth=2,
                markersize=8, label='传统混合预编码')
        plt.plot(user_range, results_uncoordinated, 'm-o', linewidth=2,
                markersize=8, label='非协调本地预编码')

        plt.xlabel('用户数 K', fontsize=14)
        plt.ylabel('和速率 (bps/Hz)', fontsize=14)
        plt.title(f'不同用户数下的和速率性能对比 (SNR = {snr_db} dB)', fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'sum_rate_vs_users.pdf'),
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir, 'sum_rate_vs_users.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"\n性能对比图已保存")

        # 保存数据
        np.savez(os.path.join(self.output_dir, 'user_sweep_data.npz'),
                user_range=user_range,
                dp_altmin=results_dp_altmin,
                fully_digital=results_fully_digital,
                traditional=results_traditional,
                uncoordinated=results_uncoordinated)

    def run_kd_sweep(self, kd_range: list = None, snr_db: float = 0,
                     num_runs: int = 10):
        """
        运行TTD组数扫描实验

        Parameters:
        -----------
        kd_range : list
            KD范围
        snr_db : float
            固定SNR (dB)
        num_runs : int
            独立运行次数
        """
        print("\n" + "=" * 80)
        print("实验4: 和速率 vs TTD组数")
        print("=" * 80)

        if kd_range is None:
            kd_range = [4, 8, 16, 32, 64]

        snr_linear = 10 ** (snr_db / 10)
        noise_power = self.config['Pmax'] / snr_linear

        results_dp_altmin = []
        results_traditional = []

        for KD in kd_range:
            print(f"\nTTD组数 KD = {KD}")

            # 更新配置
            config = self.config.copy()
            config['KD'] = KD
            system = CellFreeMIMOSystem(config)
            precoder = DistributedPrecoder(system)
            baselines = BaselineSchemes(system)

            sum_rate_dp = []
            sum_rate_th = []

            for run in range(num_runs):
                H, _ = system.generate_channel_cluster_based()

                # DP-AltMin
                res_dp = precoder.distributed_dp_altmin(H, max_iter=30)
                sum_rate_dp.append(res_dp['final_sum_rate'])

                # 传统混合（无TTD）
                res_th = baselines.traditional_hybrid(H)
                sum_rate_th.append(res_th['sum_rate'])

            results_dp_altmin.append(np.mean(sum_rate_dp))
            results_traditional.append(np.mean(sum_rate_th))

            print(f"  DP-AltMin: {results_dp_altmin[-1]:.2f} bps/Hz")
            print(f"  传统混合: {results_traditional[-1]:.2f} bps/Hz")

        # 绘制性能对比图
        plt.figure(figsize=(10, 6))
        plt.plot(kd_range, results_dp_altmin, 'b-s', linewidth=2,
                markersize=8, label='DP-AltMin (带TTD)')
        plt.plot(kd_range, results_traditional, 'r-^', linewidth=2,
                markersize=8, label='传统混合预编码 (无TTD)')

        plt.xlabel('TTD组数 KD', fontsize=14)
        plt.ylabel('和速率 (bps/Hz)', fontsize=14)
        plt.title(f'不同TTD组数下的和速率性能对比 (SNR = {snr_db} dB)', fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'sum_rate_vs_kd.pdf'),
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir, 'sum_rate_vs_kd.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"\n性能对比图已保存")

        # 保存数据
        np.savez(os.path.join(self.output_dir, 'kd_sweep_data.npz'),
                kd_range=kd_range,
                dp_altmin=results_dp_altmin,
                traditional=results_traditional)

    def run_all_experiments(self):
        """运行所有实验"""
        print("\n" + "=" * 80)
        print("开始运行所有仿真实验")
        print("=" * 80)

        start_time = time.time()

        # 实验1: 收敛性分析
        self.run_convergence_analysis(num_runs=5)

        # 实验2: SNR扫描
        self.run_snr_sweep(num_runs=5)

        # 实验3: 用户数扫描
        self.run_user_sweep(num_runs=5)

        # 实验4: TTD组数扫描
        self.run_kd_sweep(num_runs=5)

        end_time = time.time()
        total_time = end_time - start_time

        print("\n" + "=" * 80)
        print("所有实验完成!")
        print(f"总运行时间: {total_time/60:.2f} 分钟")
        print("=" * 80)

        # 保存实验配置
        with open(os.path.join(self.output_dir, 'simulation_config.json'), 'w') as f:
            json.dump(self.config, f, indent=2)

        print(f"\n结果已保存到: {self.output_dir}")


def main():
    """主函数"""
    # 创建仿真实验对象
    sim = SimulationExperiment(output_dir="../results")

    # 运行所有实验
    sim.run_all_experiments()


if __name__ == "__main__":
    main()
