"""
分布式预编码与时延-相位联合优化算法
基于WMMSE和交替优化
"""

import numpy as np
from typing import List, Tuple, Dict
import sys
sys.path.append('..')
from system_model import CellFreeMIMOSystem

class DistributedPrecoder:
    """分布式预编码器"""

    def __init__(self, system: CellFreeMIMOSystem):
        """
        初始化预编码器

        Parameters:
        -----------
        system : CellFreeMIMOSystem
            Cell-Free MIMO系统对象
        """
        self.system = system
        self.eps = 1e-10  # 数值稳定性

    def initialize_precoder(self, H: np.ndarray) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
        """
        初始化预编码矩阵

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)

        Returns:
        --------
        A : List[np.ndarray]
            相移矩阵列表
        T : List[np.ndarray]
            TTD时延矩阵列表
        D : List[np.ndarray]
            数字预编码矩阵列表
        """
        A = []
        T = []
        D = []

        for b in range(self.system.B):
            # 初始化相移矩阵A (Nt x NRF)
            # 使用单位模约束，随机相位
            A_b = np.exp(1j * 2 * np.pi * np.random.rand(self.system.Nt, self.system.NRF))
            A_b = A_b / np.abs(A_b)  # 归一化到单位模
            A.append(A_b)

            # 初始化TTD时延矩阵T (KD x NRF x M)
            # TTD结构：每个KD组处理P个天线，通过时延连接到NRF个RF链
            T_b = np.ones((self.system.KD, self.system.NRF, self.system.M), dtype=complex)
            for m in range(self.system.M):
                # 为每个TTD组生成随机时延
                t_delay = np.random.rand(self.system.KD, self.system.NRF) * 20e-9
                f = self.system.frequencies[m]
                T_b[:, :, m] = np.exp(-1j * 2 * np.pi * f * t_delay)
            T.append(T_b)

            # 初始化数字预编码矩阵D (NRF x K x M)
            D_b = np.zeros((self.system.NRF, self.system.K, self.system.M), dtype=complex)
            for m in range(self.system.M):
                # 使用随机初始化，避免数值问题
                D_b[:, :, m] = np.random.randn(self.system.NRF, self.system.K) + \
                               1j * np.random.randn(self.system.NRF, self.system.K)
                # 归一化
                D_b[:, :, m] = D_b[:, :, m] / np.linalg.norm(D_b[:, :, m], 'fro')
            D.append(D_b)

        return A, T, D

    def wmmse_update(self, H: np.ndarray, A: List[np.ndarray], T: List[np.ndarray],
                     D: List[np.ndarray], noise_power: float = 1e-10) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
        """
        WMMSE更新

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)
        A : List[np.ndarray]
            当前相移矩阵
        T : List[np.ndarray]
            当前TTD时延矩阵
        D : List[np.ndarray]
            当前数字预编码矩阵
        noise_power : float
            噪声功率

        Returns:
        --------
        u, w, D_new : 更新后的参数
        """
        u = np.zeros((self.system.K, self.system.M), dtype=complex)
        w = np.zeros((self.system.K, self.system.M))
        D_new = []

        for m in range(self.system.M):
            # 计算等效信道矩阵
            H_eff = np.zeros((self.system.K, self.system.B * self.system.NRF), dtype=complex)
            for k in range(self.system.K):
                col_idx = 0
                for b in range(self.system.B):
                    A_b = A[b]
                    T_b = T[b][:,:,m] if len(T[b].shape) == 3 else T[b]
                    D_b = D[b][:,:,m] if len(D[b].shape) == 3 else D[b]

                    # 等效信道: h^H * A * T
                    h_eff = H[b, k, m, :] @ A_b @ T_b
                    H_eff[k, col_idx:col_idx+self.system.NRF] = h_eff
                    col_idx += self.system.NRF

            # MMSE接收机
            for k in range(self.system.K):
                # 提取第k个用户的等效信道
                h_k = H_eff[k, :]

                # 计算干扰+噪声协方差矩阵
                R = noise_power * np.eye(self.system.B * self.system.NRF)
                for j in range(self.system.K):
                    if j != k:
                        h_j = H_eff[j, :]
                        # 提取第j个用户的预编码向量
                        d_j = np.zeros(self.system.B * self.system.NRF, dtype=complex)
                        col_idx = 0
                        for b in range(self.system.B):
                            D_b = D[b][:,:,m] if len(D[b].shape) == 3 else D[b]
                            d_j[col_idx:col_idx+self.system.NRF] = D_b[:, j]
                            col_idx += self.system.NRF
                        R += np.outer(h_j, h_j.conj()) * np.abs(d_j)**2

                # MMSE接收机系数
                d_k = np.zeros(self.system.B * self.system.NRF, dtype=complex)
                col_idx = 0
                for b in range(self.system.B):
                    D_b = D[b][:,:,m] if len(D[b].shape) == 3 else D[b]
                    d_k[col_idx:col_idx+self.system.NRF] = D_b[:, k]
                    col_idx += self.system.NRF

                u[k, m] = np.conj(h_k @ d_k) / (np.conj(h_k) @ R @ h_k + self.eps)

                # MSE权重
                e_k = 1 - np.real(u[k, m] * h_k @ d_k)
                w[k, m] = 1 / (e_k + self.eps)

        # 更新数字预编码矩阵
        for b in range(self.system.B):
            D_b_new = np.zeros((self.system.NRF, self.system.K, self.system.M), dtype=complex)
            for m in range(self.system.M):
                # 构建优化问题
                H_eff_b = np.zeros((self.system.K, self.system.NRF), dtype=complex)
                for k in range(self.system.K):
                    A_b = A[b]
                    T_b = T[b][:,:,m] if len(T[b].shape) == 3 else T[b]
                    H_eff_b[k, :] = H[b, k, m, :] @ A_b @ T_b

                # 加权最小二乘解
                W = np.diag(w[:, m])
                U = np.diag(u[:, m])

                # 目标函数最小化
                # min ||W^{1/2}(U H_eff D - I)||_F^2
                D_b_new[:, :, m] = np.linalg.pinv(H_eff_b) @ np.linalg.inv(U) @ np.sqrt(W)

            D_new.append(D_b_new)

        return u, w, D_new

    def update_phase_shift(self, H: np.ndarray, T: List[np.ndarray],
                           D: List[np.ndarray]) -> List[np.ndarray]:
        """
        更新相移矩阵A (在CPU端执行)

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)
        T : List[np.ndarray]
            TTD时延矩阵
        D : List[np.ndarray]
            数字预编码矩阵

        Returns:
        --------
        A_new : List[np.ndarray]
            更新后的相移矩阵
        """
        A_new = []

        for b in range(self.system.B):
            # 构建目标矩阵F
            F_b = np.zeros((self.system.Nt, self.system.NRF), dtype=complex)
            for m in range(self.system.M):
                T_b = T[b][:,:,m] if len(T[b].shape) == 3 else T[b]
                D_b = D[b][:,:,m] if len(D[b].shape) == 3 else D[b]
                F_b += H[b, :, m, :].T @ D_b @ T_b.T.conj()

            # 使用Cauchy-Schwarz不等式求解单位模约束
            # 最优解: A = (1/Nt) exp(j angle(F))
            A_b_new = np.exp(1j * np.angle(F_b)) / self.system.Nt
            A_new.append(A_b_new)

        return A_new

    def update_ttd_delay(self, H: np.ndarray, A: List[np.ndarray],
                         D: List[np.ndarray], Tmax: float = 20e-9) -> List[np.ndarray]:
        """
        更新TTD时延矩阵T (在AP端执行)

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)
        A : List[np.ndarray]
            相移矩阵
        D : List[np.ndarray]
            数字预编码矩阵
        Tmax : float
            最大时延

        Returns:
        --------
        T_new : List[np.ndarray]
            更新后的TTD时延矩阵
        """
        T_new = []

        for b in range(self.system.B):
            T_b_new = np.zeros((self.system.KD, self.system.NRF, self.system.M), dtype=complex)

            for k in range(self.system.KD):
                for n in range(self.system.NRF):
                    # 构建 A_b 的第k块
                    A_k = A[b][k*self.system.P:(k+1)*self.system.P, n]

                    # 对每个子载波构建优化目标
                    objective = np.zeros(self.system.M, dtype=complex)
                    for m in range(self.system.M):
                        # 计算V矩阵
                        V = np.zeros(self.system.Nt, dtype=complex)
                        for k_user in range(self.system.K):
                            V += H[b, k_user, m, :].conj() * D[b][n, k_user, m]

                        objective[m] = A_k @ V

                    # 一维搜索最优时延
                    # max Re{sum_m objective[m] * exp(-j2πf_m t)}
                    # 在[0, Tmax]范围内搜索
                    t_samples = np.linspace(0, Tmax, 100)
                    best_t = 0
                    best_value = -np.inf

                    for t in t_samples:
                        value = 0
                        for m in range(self.system.M):
                            f = self.system.frequencies[m]
                            value += np.real(objective[m] * np.exp(-1j * 2*np.pi * f * t))

                        if value > best_value:
                            best_value = value
                            best_t = t

                    # 设置TTD时延
                    for m in range(self.system.M):
                        f = self.system.frequencies[m]
                        T_b_new[k, n, m] = np.exp(-1j * 2*np.pi * f * best_t)

            T_new.append(T_b_new)

        return T_new

    def distributed_dp_altmin(self, H: np.ndarray, max_iter: int = 50,
                              convergence_threshold: float = 1e-4) -> Dict:
        """
        分布式DP-AltMin算法主函数

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)
        max_iter : int
            最大迭代次数
        convergence_threshold : float
            收敛阈值

        Returns:
        --------
        results : Dict
            包含优化结果和性能指标
        """
        print("\n开始分布式DP-AltMin算法...")
        print(f"最大迭代次数: {max_iter}")

        # 初始化
        A, T, D = self.initialize_precoder(H)

        # 记录收敛过程
        sum_rates = []
        iterations = []

        for iter in range(max_iter):
            # 1. CPU端：更新数字预编码矩阵D (WMMSE)
            u, w, D = self.wmmse_update(H, A, T, D)

            # 2. AP端：更新相移矩阵A
            A = self.update_phase_shift(H, T, D)

            # 3. AP端：更新TTD时延矩阵T
            T = self.update_ttd_delay(H, A, D)

            # 4. 计算当前性能
            SINR = self.system.compute_sinr(H, A, T, D)
            sum_rate = self.system.compute_sum_rate(SINR)

            sum_rates.append(sum_rate)
            iterations.append(iter + 1)

            # 打印进度
            if (iter + 1) % 5 == 0:
                print(f"迭代 {iter+1}/{max_iter}: 和速率 = {sum_rate:.2f} bps/Hz")

            # 检查收敛
            if iter > 0:
                if abs(sum_rates[-1] - sum_rates[-2]) < convergence_threshold * sum_rates[-2]:
                    print(f"\n算法在 {iter+1} 次迭代后收敛")
                    break

        print(f"\n最终和速率: {sum_rates[-1]:.2f} bps/Hz")

        return {
            'A': A,
            'T': T,
            'D': D,
            'sum_rates': sum_rates,
            'iterations': iterations,
            'final_sum_rate': sum_rates[-1]
        }


class BaselineSchemes:
    """基准方案"""

    def __init__(self, system: CellFreeMIMOSystem):
        self.system = system
        self.precoder = DistributedPrecoder(system)

    def fully_digital(self, H: np.ndarray) -> Dict:
        """
        全数字预编码基准方案

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)

        Returns:
        --------
        results : Dict
            性能结果
        """
        print("\n运行全数字预编码基准方案...")

        # 全数字预编码：每个AP使用所有天线，无RF约束
        A = []
        T = []
        D = []

        for b in range(self.system.B):
            # 全数字：A = I (单位矩阵)
            A.append(np.eye(self.system.Nt))
            # 无TTD：T = I
            T.append(np.ones((self.system.Nt, self.system.Nt, self.system.M), dtype=complex))

            # 数字预编码：使用迫零或MMSE
            D_b = np.zeros((self.system.Nt, self.system.K, self.system.M), dtype=complex)
            for m in range(self.system.M):
                H_m = H[b, :, m, :]  # K x Nt
                # MMSE预编码
                D_b[:, :, m] = H_m.conj().T @ np.linalg.inv(H_m @ H_m.conj().T +
                                     0.1 * np.eye(self.system.K))
            D.append(D_b)

        SINR = self.system.compute_sinr(H, A, T, D)
        sum_rate = self.system.compute_sum_rate(SINR)

        print(f"全数字预编码和速率: {sum_rate:.2f} bps/Hz")

        return {
            'sum_rate': sum_rate,
            'SINR': SINR
        }

    def traditional_hybrid(self, H: np.ndarray) -> Dict:
        """
        传统混合预编码基准方案（无TTD）

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵

        Returns:
        --------
        results : Dict
            性能结果
        """
        print("\n运行传统混合预编码基准方案（无TTD）...")

        # 初始化预编码器
        A, T, D = self.precoder.initialize_precoder(H)

        # 传统混合预编码：T = I (无TTD结构)
        for b in range(self.system.B):
            T[b] = np.ones((self.system.KD, self.system.NRF, self.system.M), dtype=complex)

        # 运行WMMSE优化
        for iter in range(30):
            u, w, D = self.precoder.wmmse_update(H, A, T, D)
            A = self.precoder.update_phase_shift(H, T, D)

        SINR = self.system.compute_sinr(H, A, T, D)
        sum_rate = self.system.compute_sum_rate(SINR)

        print(f"传统混合预编码和速率: {sum_rate:.2f} bps/Hz")

        return {
            'sum_rate': sum_rate,
            'SINR': SINR
        }

    def uncoordinated_local(self, H: np.ndarray) -> Dict:
        """
        非协调本地预编码基准方案

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵

        Returns:
        --------
        results : Dict
            性能结果
        """
        print("\n运行非协调本地预编码基准方案...")

        A = []
        T = []
        D = []

        # 每个AP独立优化，不考虑其他AP的干扰
        for b in range(self.system.B):
            A_b = np.exp(1j * 2 * np.pi * np.random.rand(self.system.Nt, self.system.NRF))
            A_b = A_b / np.abs(A_b)
            A.append(A_b)

            T_b = np.ones((self.system.KD, self.system.NRF, self.system.M), dtype=complex)
            T.append(T_b)

            D_b = np.zeros((self.system.NRF, self.system.K, self.system.M), dtype=complex)
            for m in range(self.system.M):
                H_m = H[b, :, m, :] @ A_b  # K x NRF
                D_b[:, :, m] = H_m.conj().T @ np.linalg.inv(H_m @ H_m.conj().T +
                                     0.1 * np.eye(self.system.K))
            D.append(D_b)

        SINR = self.system.compute_sinr(H, A, T, D)
        sum_rate = self.system.compute_sum_rate(SINR)

        print(f"非协调本地预编码和速率: {sum_rate:.2f} bps/Hz")

        return {
            'sum_rate': sum_rate,
            'SINR': SINR
        }


if __name__ == "__main__":
    # 测试预编码算法
    from system_model import CellFreeMIMOSystem

    config = {
        'B': 4,
        'K': 4,
        'Nt': 256,
        'NRF': 4,
        'KD': 16,
        'M': 128,
        'fc': 300e9,
        'W': 30e9
    }

    system = CellFreeMIMOSystem(config)
    print("\n生成信道...")
    H, _ = system.generate_channel_cluster_based()

    print("\n测试分布式DP-AltMin算法...")
    precoder = DistributedPrecoder(system)
    results = precoder.distributed_dp_altmin(H, max_iter=20)

    print(f"\n最终和速率: {results['final_sum_rate']:.2f} bps/Hz")
