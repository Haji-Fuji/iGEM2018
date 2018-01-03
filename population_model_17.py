import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

np.set_printoptions(threshold=np.inf)


def func(v, t, a, b, c, e, f, K1, K2, ki, kii, u, d1, d2, d3, d4, d5):
    coli = v[0]
    C8 = v[1]
    IP = v[2]
    MazF = v[3]
    DiMazF = v[4] #MazFはホモダイマーになってmRNAを分解します

    return [((a *K1 * coli) / (K1 + DiMazF)) - f * coli,
            b * coli - d2 * C8 - f * C8,
            c * u * C8 - d3 * IP - f * IP,
            e * IP/(K2+MazF) - 2 * ki * MazF + 2 * kii * DiMazF - d4 * MazF,
            ki * MazF - kii * DiMazF - d5 * DiMazF]


# 初期値
initial = [0.1, 0, 0, 0, 0]

# パラメータ
a = 0.0221  # 大腸菌増殖速度(実験値)
b = 0.001 # C8合成速度(実験値)
c = 0.02    # iP合成速度
e = 0.03    # MazF合成速度
f = 0.001   # 流失速度(varied)
K1 = 0.46   # DiMazFの増殖速度への影響(去年のデータ)
K2 = 1 # MazFのミカエリス定数
ki = 6.82 #Formation rate of MazF dimer(去年のデータ)
kii = 6.24 #Dissociation rate of MazF dimer(去年のデータ)
u = 2       # ヒト細胞の数
d1 = 0.46   # MazFによる抑制
d2 = 0.00088   # C8の分解速度(去年の参考文献より)
d3 = 0.06   # iPの分解速度(#Cytokinin oxidase and the regulation of cytokinin degradation)
d4 = 0.7   # MazFの分解速度(去年のデータより)
d5 = 0.17 # MazF dimerの分解速度(去年のデータ)

# 時間
t = np.arange(0, 1440 * 5, 0.1, dtype=float)  # 1440分＝1day

# 時間測定開始
start_time = time.time()

# 計算
v0 = initial
v = odeint(func, v0, t, args=(a, b, c, e, f, K1, K2, ki, kii, u, d1, d2, d3, d4, d5))

# 表示
plt.plot(t, v[:, 0], label='E. coli')
plt.plot(t, v[:, 1], label='C8')
plt.plot(t, v[:, 2], label='IP')
plt.plot(t, v[:, 3], label='MazF')
plt.plot(t, v[:, 4], label='DiMazF')
plt.xlabel('t [min]')  # X軸
plt.legend()
plt.show()
# plt.savefig('aaa.png')
plt.close()

# 時間測定終了
end_time = time.time()
print(str(end_time - start_time) + '秒')