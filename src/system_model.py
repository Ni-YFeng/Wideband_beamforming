"""
Cell-Free MIMO系统模型定义
包含系统参数、信道模型和预编码结构
"""

import numpy as np
from typing import List, Tuple, Dict
import json

class CellFreeMIMOSystem:
    """Cell-Free MIMO系统类"""

    def __init__(self, config: Dict):
        """
        初始化系统参数

        Parameters:
        -----------
        config : Dict
            系统配置字典，包含：
            - B: AP数量
            - K: 用户数量
            - Nt: 每个AP的天线数
            - NRF: RF链数
            - KD: TTD组数
            - M: OFDM子载波数
            - fc: 载波频率 (Hz)
            - W: 带宽 (Hz)
            - Pmax: 最大发射功率
        """
        self.B = config.get('B', 4)  # AP数量
        self.K = config.get('K', 4)  # 用户数量
        self.Nt = config.get('Nt', 256)  # 每个AP的天线数
        self.NRF = config.get('NRF', 4)  # RF链数
        self.KD = config.get('KD', 16)  # TTD组数
        self.M = config.get('M', 128)  # OFDM子载波数
        self.fc = config.get('fc', 300e9)  # 载波频率 300 GHz
        self.W = config.get('W', 30e9)  # 带宽 30 GHz
        self.Pmax = config.get('Pmax', 1.0)  # 归一化功率

        # 计算其他参数
        self.P = self.Nt // self.KD  # 每个TTD组的天线数
        self.c = 3e8  # 光速

        # 频率数组
        self.frequencies = self._generate_frequencies()

        print(f"系统初始化完成:")
        print(f"  AP数量 B = {self.B}")
        print(f"  用户数量 K = {self.K}")
        print(f"  天线数 Nt = {self.Nt}")
        print(f"  RF链数 NRF = {self.NRF}")
        print(f"  TTD组数 KD = {self.KD}")
        print(f"  子载波数 M = {self.M}")
        print(f"  载波频率 fc = {self.fc/1e9:.1f} GHz")
        print(f"  带宽 W = {self.W/1e9:.1f} GHz")

    def _generate_frequencies(self) -> np.ndarray:
        """生成OFDM子载波频率"""
        frequencies = np.zeros(self.M)
        for m in range(self.M):
            frequencies[m] = self.fc - self.W/2 + (m) * self.W/self.M
        return frequencies

    def array_response(self, theta: float, N: int, lambda_c: float, lambda_f: float) -> np.ndarray:
        """
        计算阵列响应向量

        Parameters:
        -----------
        theta : float
            角度（弧度）
        N : int
            天线数
        lambda_c : float
            载波波长
        lambda_f : float
            当前频率的波长

        Returns:
        --------
        np.ndarray
            阵列响应向量
        """
        a = np.zeros(N, dtype=complex)
        for n in range(N):
            a[n] = np.exp(1j * 2 * np.pi * n * lambda_c/lambda_f * np.sin(theta))
        a = a / np.sqrt(N)
        return a

    def generate_channel_cluster_based(self,
                                       Lc: int = 5,
                                       Lp: int = 20,
                                       sigma_bs: float = 0.1,
                                       sigma_mt: float = 0.1,
                                       sigma_t: float = 5e-9,
                                       Tmax: float = 20e-9) -> Tuple[np.ndarray, np.ndarray]:
        """
        基于簇的宽带信道生成

        Parameters:
        -----------
        Lc : int
            簇的数量
        Lp : int
            每个簇中的路径数
        sigma_bs : float
            基站端的角度扩展
        sigma_mt : float
            移动端的角度扩展
        sigma_t : float
            时延扩展
        Tmax : float
            最大时延

        Returns:
        --------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)
        channel_params : np.ndarray
            信道参数
        """
        H = np.zeros((self.B, self.K, self.M, self.Nt), dtype=complex)

        # 存储信道参数
        channel_params = {
            'theta': [],  # 离开角
            'alpha': [],  # 到达角
            'beta': [],   # 复增益
            'delay': [],  # 时延
            'Lc': Lc,
            'Lp': Lp
        }

        for b in range(self.B):
            for k in range(self.K):
                # 为每个AP-用户对生成信道
                # 离开角 (AOD)
                Theta_c = np.random.uniform(-np.pi/2, np.pi/2, Lc)
                sigma = sigma_bs
                dTheta = -sigma/np.sqrt(2) * np.sign(np.random.randn(Lc, Lp) - 0.5) * \
                         np.log(1 - 2*np.abs(np.random.randn(Lc, Lp) - 0.5))

                # 到达角 (AOA)
                Alpha_c = np.random.uniform(-np.pi/2, np.pi/2, Lc)
                sigma = sigma_mt
                dAlpha = -sigma/np.sqrt(2) * np.sign(np.random.randn(Lc, Lp) - 0.5) * \
                         np.log(1 - 2*np.abs(np.random.randn(Lc, Lp) - 0.5))

                # 复增益
                Beta = (np.random.randn(Lc, Lp) + 1j*np.random.randn(Lc, Lp)) / np.sqrt(2)

                # 时延
                Delay_c = np.random.uniform(0, Tmax, Lc)
                sigma = sigma_t
                dDelay = -sigma/np.sqrt(2) * np.sign(np.random.randn(Lc, Lp) - 0.5) * \
                         np.log(1 - 2*np.abs(np.random.randn(Lc, Lp) - 0.5))

                # 对每个子载波计算信道
                for m in range(self.M):
                    f = self.frequencies[m]
                    lambda_f = self.c / f
                    lambda_c = self.c / self.fc

                    h = np.zeros(self.Nt, dtype=complex)
                    for c in range(Lc):
                        for p in range(Lp):
                            theta = Theta_c[c] + dTheta[c, p]
                            alpha = Alpha_c[c] + dAlpha[c, p]
                            delay = Delay_c[c] + dDelay[c, p]
                            beta = Beta[c, p]

                            # 阵列响应
                            a = self.array_response(theta, self.Nt, lambda_c, lambda_f)

                            # 信道
                            h += beta * np.exp(-1j * 2*np.pi * delay * f) * a

                    H[b, k, m, :] = h * np.sqrt(self.Nt / (Lc * Lp))

                channel_params['theta'].append(Theta_c)
                channel_params['alpha'].append(Alpha_c)
                channel_params['beta'].append(Beta)
                channel_params['delay'].append(Delay_c)

        return H, channel_params

    def compute_sinr(self,
                     H: np.ndarray,
                     A: List[np.ndarray],
                     T: List[np.ndarray],
                     D: List[np.ndarray],
                     noise_power: float = 1e-10) -> np.ndarray:
        """
        计算SINR

        Parameters:
        -----------
        H : np.ndarray
            信道矩阵 (B, K, M, Nt)
        A : List[np.ndarray]
            相移矩阵列表 (B个, 每个形状为 Nt x NRF)
        T : List[np.ndarray]
            TTD时延矩阵列表 (B个, 每个形状为 KD x NRF x M)
        D : List[np.ndarray]
            数字预编码矩阵列表 (B个, 每个形状为 NRF x K x M)
        noise_power : float
            噪声功率

        Returns:
        --------
        SINR : np.ndarray
            SINR矩阵 (K, M)
        """
        SINR = np.zeros((self.K, self.M))

        for m in range(self.M):
            for k in range(self.K):
                # 计算等效信道
                h_eff = np.zeros(self.NRF, dtype=complex)
                for b in range(self.B):
                    # TTD预编码: A * T[:,:,m] * D[:,:,m]
                    A_b = A[b]
                    T_b = T[b][:,:,m] if len(T[b].shape) == 3 else T[b]
                    D_b = D[b][:,:,m] if len(D[b].shape) == 3 else D[b]

                    F_b = A_b @ T_b @ D_b
                    h_eff += H[b, k, m, :] @ F_b

                # 信号功率
                signal_power = np.abs(h_eff[k])**2

                # 干扰功率
                interference_power = 0
                for j in range(self.K):
                    if j != k:
                        h_j = np.zeros(self.NRF, dtype=complex)
                        for b in range(self.B):
                            A_b = A[b]
                            T_b = T[b][:,:,m] if len(T[b].shape) == 3 else T[b]
                            D_b = D[b][:,:,m] if len(D[b].shape) == 3 else D[b]
                            F_b = A_b @ T_b @ D_b
                            h_j += H[b, k, m, :] @ F_b
                        interference_power += np.abs(h_j[j])**2

                SINR[k, m] = signal_power / (interference_power + noise_power)

        return SINR

    def compute_sum_rate(self, SINR: np.ndarray) -> float:
        """
        计算和速率

        Parameters:
        -----------
        SINR : np.ndarray
            SINR矩阵 (K, M)

        Returns:
        --------
        float
            和速率 (bps/Hz)
        """
        sum_rate = 0
        for m in range(self.M):
            for k in range(self.K):
                sum_rate += np.log2(1 + SINR[k, m])
        return sum_rate

    def save_config(self, filepath: str):
        """保存系统配置"""
        config = {
            'B': self.B,
            'K': self.K,
            'Nt': self.Nt,
            'NRF': self.NRF,
            'KD': self.KD,
            'M': self.M,
            'fc': self.fc,
            'W': self.W,
            'Pmax': self.Pmax
        }
        with open(filepath, 'w') as f:
            json.dump(config, f, indent=2)
        print(f"配置已保存到 {filepath}")


if __name__ == "__main__":
    # 测试系统初始化
    config = {
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

    system = CellFreeMIMOSystem(config)
    print("\n生成测试信道...")
    H, params = system.generate_channel_cluster_based()
    print(f"信道形状: {H.shape}")
    print(f"信道功率: {np.mean(np.abs(H)):.6f}")
