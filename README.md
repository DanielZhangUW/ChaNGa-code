Single-Particle Hill Test (ChaNGa vs Quinn)
===========================================

Purpose
-------
This repo contains a Python reference integrator for Hill's equations and a
ChaNGa sliding-patch run card. The main goal is to compare a single-particle
orbit step-by-step in code units (Omega=1, 2*pi code time per orbit).

Key Files
---------
- third_party/ChaNGa/single_particle_hill.param
  ChaNGa run card (sliding patch, 100 steps per orbit).
- tools/make_single_particle_tipsy.py
  Generates third_party/ChaNGa/single_particle.std from a chosen IC.
- task2/plot_task2.py
  Quinn/Hill reference integrator (Python).
- task2/plot_time_series.py
  Plots x(t), y(t), vx(t), vy(t) for ChaNGa vs Quinn.

Quick Start (Epicycle IC)
-------------------------
From the repo root:

  # Defaults are epicycle (x0=1, y0=0, vx0=0, vy0=-2).
  python3 tools/make_single_particle_tipsy.py --add-guard
  cd third_party/ChaNGa
  ./ChaNGa single_particle_hill.param +p1
  cd ..
  python3 task2/plot_time_series.py --output task2/single_particle_timeseries.png

Expected Outputs
----------------
- third_party/ChaNGa/single_particle.std
  Tipsy input with one physical particle + guard particle.
- third_party/ChaNGa/single_particle_hill.000001 ... .000100
  ChaNGa snapshots (one per step).
- task2/single_particle_timeseries.png
  4-panel plot comparing Quinn vs ChaNGa time series.

  std output:

python3 - <<'PY'
import pynbody as pb
s = pb.load("third_party/ChaNGa/single_particle.std")
p = s[s["mass"]==s["mass"].max()][0]
print("mass =", p["mass"])
print("pos  =", p["pos"])
print("vel  =", p["vel"])
PY

mass = [1.]
pos  = [[1. 0. 0.]]
vel  = [[ 0. -2.  0.]]


Units and Time Step
-------------------
- Code units are chosen so Omega = 1 and one orbit is 2*pi in code time.
- Length unit: dKpcUnit = 4.848e-9 so 1 code length ~ 1 AU.
- Time step: dDelta = 2*pi/100 (100 steps per orbit).
- dMsolUnit = 1.0, one code mass unit = 1 Msun

Consistency Checks
------------------
- The IC defaults are consistent across:
  - tools/make_single_particle_tipsy.py
  - task2/plot_task2.py
  - task2/plot_time_series.py
- plot_time_series.py prints:
  - Quinn requested IC (code units)
  - ChaNGa first snapshot (code units)
  This is useful to confirm the run used the expected IC.

Notes
-----
- Units: dKpcUnit = 4.848e-9 so 1 code length unit ~ 1 AU.
- ChaNGa reads third_party/ChaNGa/single_particle.std as input.
- Plot script prints both Quinn requested IC and ChaNGa first snapshot IC.

1) README.md 

README

就是说明书：告诉你 Task2 的目标（Quinn(Python) vs ChaNGa 对 Hill’s equations 的对比）和运行流程。

2) tools/make_single_particle_tipsy.py

它负责生成 ChaNGa 要读的初始条件文件（Tipsy 格式）。

输入：--x0 --y0 --vx0 --vy0（这些是 code units）

输出：third_party/ChaNGa/single_particle.std

还会加一个“guard particle”（质量10e−12
）来避免 ChaNGa 的边界/奇异情况。

3) single_particle.std

就是第 2 步生成的 ChaNGa 初始条件二进制文件（你不用读懂内容，知道它是“IC 文件”就够）。

4) single_particle_hill.param

ChaNGa 的参数文件：控制

读哪个 IC (achInFile)

输出叫什么 (achOutName)

多少步 (nSteps) / 步长 dt (dDelta)

以及最重要的：patch/periodic 这些设置（会导致 y 坐标“折回/跳变”）

5) plot_task2.py

这是你的 Python 版 Quinn integrator（参考答案）。
最关键函数：run_quinn_states_for_ic() 会给出每一步的 x,y,vx,vy。

这里面有一堆你问的变量：

vx, vy 就是最朴素的速度：dx/dt, dy/dt

px = vx + y, py = vy - x（类似 canonical momentum 的组合）

Px = vx - y, Py = vy + 2x（特别重要：Py=vy+2x 在“无外力 Hill 方程”里是守恒量）

所以：如果初始是 (1,0,0,-2)，那么

Py(0)= -2 + 2*1 = 0

并且理论上 Py(t) 永远等于 0
这就是你后来看到“蓝线一直是 0”为什么其实可能是对的（它在检验守恒量）。

6) plot_time_series.py 

plot_time_series

它做的事是：

读 ChaNGa 输出的 snapshots

找到真正的物理粒子（按质量筛掉 guard）

把 pos/vel 转回 code units

再用 Python Quinn 算一条曲线叠上去

画 4 张：x(t), y(t), vx(t), vy(t)

教授很可能会问：y(t) 为啥有竖直尖峰？
答案：因为 patch/periodic 条件下，ChaNGa 可能把 y “wrap”到一个区间里，过边界就会瞬间跳变。这不是物理突然跳，是坐标表示方式的问题。要么“unwrap y”，要么把 box 调大/振幅调小，让它不跨边界。

7) compare_single_orbit.py

画的是 x–y 平面轨道散点图（100 steps）：

蓝：Python Quinn

红：ChaNGa

如果 ChaNGa 的红线变成一条斜线，很常见就是 y 被 wrap 了，导致点被“折回”。

Result of epicycle:
ChaNGa INPUT single_particle.std (raw units):
  mass = 1.000000e+00
  pos  = [[1. 0. 0.]]
  vel  = [[ 0. -2.  0.]]
ChaNGa OUTPUT single_particle_hill.000001:
  mass (Msol) = 1.000000e+00
  pos (kpc) = [[ 4.83843358e-09 -9.56642186e-12  0.00000000e+00]]
  vel (km/s) = [[-1.87046659 -1.8558728   0.        ]]
Quinn IC (code units): x0=1.000000, y0=0.000000, vx0=0.000000, vy0=-2.000000
Per-step diagnostics (kpc, radians):
Step 0, t=0.000000:
  Quinn : x=4.838430e-09, y=-6.086164e-10, r=4.876559e-09, phi=-1.251308e-01
  ChaNGa: x=4.838434e-09, y=-9.566422e-12, r=4.838443e-09, phi=-1.977171e-03
  delta : dx=3.148826e-15, dy=5.990500e-10, dr=-3.811548e-11, raw_dphi=1.231536e-01, wrapped_dphi=1.231536e-01 (deg=7.056)
Step 1, t=0.062832:
  Quinn : x=4.809760e-09, y=-1.214830e-09, r=4.960806e-09, phi=-2.474017e-01
  ChaNGa: x=4.809772e-09, y=-3.822793e-11, r=4.809924e-09, phi=-7.947804e-03
  delta : dx=1.267036e-14, dy=1.176602e-09, dr=-1.508821e-10, raw_dphi=2.394539e-01, wrapped_dphi=2.394539e-01 (deg=13.720)
Step 2, t=0.125664:
  Quinn : x=4.762100e-09, y=-1.816248e-09, r=5.096701e-09, phi=-3.643666e-01
  ChaNGa: x=4.762128e-09, y=-8.587142e-11, r=4.762903e-09, phi=-1.803020e-02
  delta : dx=2.807201e-14, dy=1.730376e-09, dr=-3.337979e-10, raw_dphi=3.463364e-01, wrapped_dphi=3.463364e-01 (deg=19.844)
Step 3, t=0.188496:
  Quinn : x=4.695641e-09, y=-2.410495e-09, r=5.278213e-09, phi=-4.742684e-01
  ChaNGa: x=4.695691e-09, y=-1.523088e-10, r=4.698161e-09, phi=-3.242451e-02
  delta : dx=4.989625e-14, dy=2.258186e-09, dr=-5.800525e-10, raw_dphi=4.418439e-01, wrapped_dphi=4.418439e-01 (deg=25.316)
Step 4, t=0.251327:
  Quinn : x=4.610645e-09, y=-2.995226e-09, r=5.498129e-09, phi=-5.761171e-01
  ChaNGa: x=4.610722e-09, y=-2.372780e-10, r=4.616823e-09, phi=-5.141687e-02
  delta : dx=7.757380e-14, dy=2.757948e-09, dr=-8.813056e-10, raw_dphi=5.247002e-01, wrapped_dphi=5.247002e-01 (deg=30.063)
Step 5, t=0.314159:
  Quinn : x=4.507446e-09, y=-3.568133e-09, r=5.748795e-09, phi=-6.696033e-01
  ChaNGa: x=4.507556e-09, y=-3.404436e-10, r=4.520395e-09, phi=-7.538419e-02
  delta : dx=1.107483e-13, dy=3.227689e-09, dr=-1.228400e-09, raw_dphi=5.942191e-01, wrapped_dphi=5.942191e-01 (deg=34.046)
Step 6, t=0.376991:
  Quinn : x=4.386452e-09, y=-4.126953e-09, r=6.022682e-09, phi=-7.549264e-01
  ChaNGa: x=4.386602e-09, y=-4.613985e-10, r=4.410801e-09, phi=-1.047982e-01
  delta : dx=1.494453e-13, dy=3.665554e-09, dr=-1.611882e-09, raw_dphi=6.501281e-01, wrapped_dphi=6.501281e-01 (deg=37.250)
Step 7, t=0.439823:
  Quinn : x=4.248142e-09, y=-4.669480e-09, r=6.312745e-09, phi=-8.326109e-01
  ChaNGa: x=4.248335e-09, y=-5.996653e-10, r=4.290448e-09, phi=-1.402266e-01
  delta : dx=1.930799e-13, dy=4.069815e-09, dr=-2.022297e-09, raw_dphi=6.923843e-01, wrapped_dphi=6.923843e-01 (deg=39.671)
Step 8, t=0.502655:
  Quinn : x=4.093060e-09, y=-5.193574e-09, r=6.612590e-09, phi=-9.033532e-01
  ChaNGa: x=4.093302e-09, y=-7.546983e-10, r=4.162294e-09, phi=-1.823264e-01
  delta : dx=2.417645e-13, dy=4.438875e-09, dr=-2.450296e-09, raw_dphi=7.210268e-01, wrapped_dphi=7.210268e-01 (deg=41.312)
Step 9, t=0.565487:
  Quinn : x=3.921820e-09, y=-5.697163e-09, r=6.916527e-09, phi=-9.679109e-01
  ChaNGa: x=3.922114e-09, y=-9.258857e-10, r=4.029919e-09, phi=-2.318238e-01
  delta : dx=2.947082e-13, dy=4.771278e-09, dr=-2.886608e-09, raw_dphi=7.360871e-01, wrapped_dphi=7.360871e-01 (deg=42.175)
Step 10, t=0.628319:
  Quinn : x=3.735097e-09, y=-6.178262e-09, r=7.219547e-09, phi=-1.027035e+00
  ChaNGa: x=3.735448e-09, y=-1.112552e-09, r=3.897608e-09, phi=-2.894705e-01
  delta : dx=3.514774e-13, dy=5.065710e-09, dr=-3.321940e-09, raw_dphi=7.375640e-01, wrapped_dphi=7.375640e-01 (deg=42.259)
Step 11, t=0.691150:
  Quinn : x=3.533628e-09, y=-6.634969e-09, r=7.517270e-09, phi=-1.081428e+00
  ChaNGa: x=3.534040e-09, y=-1.313960e-09, r=3.770401e-09, phi=-3.559633e-01
  delta : dx=4.114868e-13, dy=5.321009e-09, dr=-3.746868e-09, raw_dphi=7.254651e-01, wrapped_dphi=7.254651e-01 (deg=41.566)
Step 12, t=0.753982:
  Quinn : x=3.318209e-09, y=-7.065483e-09, r=7.805867e-09, phi=-1.131733e+00
  ChaNGa: x=3.318684e-09, y=-1.529316e-09, r=3.654103e-09, phi=-4.318153e-01
  delta : dx=4.745973e-13, dy=5.536167e-09, dr=-4.151764e-09, raw_dphi=6.999178e-01, wrapped_dphi=6.999178e-01 (deg=40.102)
Step 13, t=0.816814:
  Quinn : x=3.089691e-09, y=-7.468103e-09, r=8.082002e-09, phi=-1.178520e+00
  ChaNGa: x=3.090231e-09, y=-1.757769e-09, r=3.555176e-09, phi=-5.171735e-01
  delta : dx=5.402106e-13, dy=5.710334e-09, dr=-4.526825e-09, raw_dphi=6.613466e-01, wrapped_dphi=6.613466e-01 (deg=37.892)
Step 14, t=0.879646:
  Quinn : x=2.848975e-09, y=-7.841240e-09, r=8.342764e-09, phi=-1.222294e+00
  ChaNGa: x=2.849582e-09, y=-1.998418e-09, r=3.480487e-09, phi=-6.115992e-01
  delta : dx=6.077371e-13, dy=5.842823e-09, dr=-4.862276e-09, raw_dphi=6.106948e-01, wrapped_dphi=6.106948e-01 (deg=34.990)
Step 15, t=0.942478:
  Quinn : x=2.597011e-09, y=-8.183422e-09, r=8.585619e-09, phi=-1.263499e+00
  ChaNGa: x=2.597687e-09, y=-2.250312e-09, r=3.436842e-09, phi=-7.138671e-01
  delta : dx=6.764010e-13, dy=5.933109e-09, dr=-5.148777e-09, raw_dphi=5.496317e-01, wrapped_dphi=5.496317e-01 (deg=31.492)
Step 16, t=1.005310:
  Quinn : x=2.334795e-09, y=-8.493296e-09, r=8.808368e-09, phi=-1.302524e+00
  ChaNGa: x=2.335541e-09, y=-2.512459e-09, r=3.430335e-09, phi=-8.218750e-01
  delta : dx=7.459575e-13, dy=5.980837e-09, dr=-5.378033e-09, raw_dphi=4.806495e-01, wrapped_dphi=4.806495e-01 (deg=27.539)
Step 17, t=1.068142:
  Quinn : x=2.063362e-09, y=-8.769640e-09, r=9.009110e-09, phi=-1.339715e+00
  ChaNGa: x=2.064177e-09, y=-2.783823e-09, r=3.465616e-09, phi=-9.327638e-01
  delta : dx=8.155633e-13, dy=5.985818e-09, dr=-5.543493e-09, raw_dphi=4.069508e-01, wrapped_dphi=4.069508e-01 (deg=23.317)
Step 18, t=1.130973:
  Quinn : x=1.783782e-09, y=-9.011364e-09, r=9.186215e-09, phi=-1.375375e+00
  ChaNGa: x=1.784667e-09, y=-3.063333e-09, r=3.545285e-09, phi=-1.043277e+00
  delta : dx=8.845401e-13, dy=5.948030e-09, dr=-5.640930e-09, raw_dphi=3.320977e-01, wrapped_dphi=3.320977e-01 (deg=19.028)
Step 19, t=1.193805:
  Quinn : x=1.497161e-09, y=-9.217511e-09, r=9.338309e-09, phi=-1.409777e+00
  ChaNGa: x=1.498113e-09, y=-3.349887e-09, r=3.669616e-09, phi=-1.150262e+00
  delta : dx=9.524150e-13, dy=5.867625e-09, dr=-5.668692e-09, raw_dphi=2.595144e-01, wrapped_dphi=2.595144e-01 (deg=14.869)
Step 20, t=1.256637:
  Quinn : x=1.204629e-09, y=-9.387270e-09, r=9.464247e-09, phi=-1.443168e+00
  ChaNGa: x=1.205647e-09, y=-3.642353e-09, r=3.836707e-09, phi=-1.251140e+00
  delta : dx=1.018366e-12, dy=5.744917e-09, dr=-5.627540e-09, raw_dphi=1.920279e-01, wrapped_dphi=1.920279e-01 (deg=11.002)
Step 21, t=1.319469:
  Quinn : x=9.073413e-10, y=-9.519969e-09, r=9.563110e-09, phi=-1.475774e+00
  ChaNGa: x=9.084232e-10, y=-3.939577e-09, r=4.042957e-09, phi=-1.344169e+00
  delta : dx=1.081858e-12, dy=5.580392e-09, dr=-5.520154e-09, raw_dphi=1.316055e-01, wrapped_dphi=1.316055e-01 (deg=7.540)
Step 22, t=1.382301:
  Quinn : x=6.064716e-10, y=-9.615085e-09, r=9.634192e-09, phi=-1.507805e+00
  ChaNGa: x=6.076134e-10, y=-4.240387e-09, r=4.283699e-09, phi=-1.428473e+00
  delta : dx=1.141717e-12, dy=5.374698e-09, dr=-5.350494e-09, raw_dphi=7.933159e-02, wrapped_dphi=7.933159e-02 (deg=4.545)
Step 23, t=1.445133:
  Quinn : x=3.032077e-10, y=-9.672241e-09, r=9.676993e-09, phi=-1.539458e+00
  ChaNGa: x=3.044056e-10, y=-4.543594e-09, r=4.553780e-09, phi=-1.503900e+00
  delta : dx=1.197904e-12, dy=5.128647e-09, dr=-5.123213e-09, raw_dphi=3.555871e-02, wrapped_dphi=3.555871e-02 (deg=2.037)
Step 24, t=1.507964:
  Quinn : x=-1.253210e-12, y=-9.691214e-09, r=9.691214e-09, phi=-1.570926e+00
  ChaNGa: x=-5.246832e-15, y=-4.848005e-09, r=4.848005e-09, phi=-1.570797e+00
  delta : dx=1.247963e-12, dy=4.843209e-09, dr=-4.843209e-09, raw_dphi=1.282318e-04, wrapped_dphi=1.282318e-04 (deg=0.007)
Step 25, t=1.570796:
  Quinn : x=-3.057092e-10, y=-9.671927e-09, r=9.676757e-09, phi=-1.602394e+00
  ChaNGa: x=-3.046139e-10, y=-5.152614e-09, r=5.161610e-09, phi=-1.629846e+00
  delta : dx=1.095262e-12, dy=4.519313e-09, dr=-4.515146e-09, raw_dphi=-2.745223e-02, wrapped_dphi=-2.745223e-02 (deg=-1.573)
Step 26, t=1.633628:
  Quinn : x=-6.089583e-10, y=-9.614456e-09, r=9.633722e-09, phi=-1.634050e+00
  ChaNGa: x=-6.078210e-10, y=-5.455821e-09, r=5.489575e-09, phi=-1.681747e+00
  delta : dx=1.137240e-12, dy=4.158635e-09, dr=-4.144148e-09, raw_dphi=-4.769700e-02, wrapped_dphi=-4.769700e-02 (deg=-2.733)
Step 27, t=1.696460:
  Quinn : x=-9.098033e-10, y=-9.519030e-09, r=9.562409e-09, phi=-1.666084e+00
  ChaNGa: x=-9.086293e-10, y=-5.756630e-09, r=5.827898e-09, phi=-1.727345e+00
  delta : dx=1.173977e-12, dy=3.762400e-09, dr=-3.734512e-09, raw_dphi=-6.126110e-02, wrapped_dphi=-6.126110e-02 (deg=-3.510)
Step 28, t=1.759292:
  Quinn : x=-1.207057e-09, y=-9.386024e-09, r=9.463320e-09, phi=-1.698696e+00
  ChaNGa: x=-1.205851e-09, y=-6.053851e-09, r=6.172778e-09, phi=-1.767410e+00
  delta : dx=1.205349e-12, dy=3.332172e-09, dr=-3.290541e-09, raw_dphi=-6.871465e-02, wrapped_dphi=-6.871465e-02 (deg=-3.937)
Step 29, t=1.822124:
  Quinn : x=-1.499545e-09, y=-9.215963e-09, r=9.337163e-09, phi=-1.732094e+00
  ChaNGa: x=-1.498314e-09, y=-6.346314e-09, r=6.520786e-09, phi=-1.802643e+00
  delta : dx=1.230463e-12, dy=2.869649e-09, dr=-2.816377e-09, raw_dphi=-7.054847e-02, wrapped_dphi=-7.054847e-02 (deg=-4.042)
Step 30, t=1.884956:
  Quinn : x=-1.786113e-09, y=-9.009519e-09, r=9.184859e-09, phi=-1.766506e+00
  ChaNGa: x=-1.784864e-09, y=-6.632864e-09, r=6.868815e-09, phi=-1.833664e+00
  delta : dx=1.248885e-12, dy=2.376655e-09, dr=-2.316043e-09, raw_dphi=-6.715749e-02, wrapped_dphi=-6.715749e-02 (deg=-3.848)
Step 31, t=1.947787:
  Quinn : x=-2.065629e-09, y=-8.767507e-09, r=9.007552e-09, phi=-1.802177e+00
  ChaNGa: x=-2.064369e-09, y=-6.912370e-09, r=7.214047e-09, phi=-1.861013e+00
  delta : dx=1.259931e-12, dy=1.855137e-09, dr=-1.793505e-09, raw_dphi=-5.883551e-02, wrapped_dphi=-5.883551e-02 (deg=-3.371)
Step 32, t=2.010619:
  Quinn : x=-2.336991e-09, y=-8.490882e-09, r=8.806623e-09, phi=-1.839381e+00
  ChaNGa: x=-2.335728e-09, y=-7.183728e-09, r=7.553911e-09, phi=-1.885156e+00
  delta : dx=1.263376e-12, dy=1.307154e-09, dr=-1.252712e-09, raw_dphi=-4.577480e-02, wrapped_dphi=-4.577480e-02 (deg=-2.623)
Step 33, t=2.073451:
  Quinn : x=-2.599127e-09, y=-8.180737e-09, r=8.583700e-09, phi=-1.878423e+00
  ChaNGa: x=-2.597868e-09, y=-7.445869e-09, r=7.886056e-09, phi=-1.906491e+00
  delta : dx=1.258824e-12, dy=7.348679e-10, dr=-6.976444e-10, raw_dphi=-2.806800e-02, wrapped_dphi=-2.806800e-02 (deg=-1.608)
Step 34, t=2.136283:
  Quinn : x=-2.851002e-09, y=-7.838295e-09, r=8.340688e-09, phi=-1.919648e+00
  ChaNGa: x=-2.849756e-09, y=-7.697756e-09, r=8.208323e-09, phi=-1.925358e+00
  delta : dx=1.245891e-12, dy=1.405386e-10, dr=-1.323654e-10, raw_dphi=-5.709827e-03, wrapped_dphi=-5.709827e-03 (deg=-0.327)
Step 35, t=2.199115:
  Quinn : x=-3.091622e-09, y=-7.464909e-09, r=8.079789e-09, phi=-1.963445e+00
  ChaNGa: x=-3.090397e-09, y=-7.938397e-09, r=8.518727e-09, phi=-1.942042e+00
  delta : dx=1.224395e-12, dy=-4.734886e-10, dr=4.389381e-10, raw_dphi=2.140216e-02, wrapped_dphi=2.140216e-02 (deg=1.226)
Step 36, t=2.261947:
  Quinn : x=-3.320036e-09, y=-7.062052e-09, r=7.803539e-09, phi=-2.010258e+00
  ChaNGa: x=-3.318842e-09, y=-8.166842e-09, r=8.815442e-09, phi=-1.956791e+00
  delta : dx=1.193940e-12, dy=-1.104790e-09, dr=1.011903e-09, raw_dphi=5.346761e-02, wrapped_dphi=5.346761e-02 (deg=3.063)
Step 37, t=2.324779:
  Quinn : x=-3.535344e-09, y=-6.631316e-09, r=7.514853e-09, phi=-2.060594e+00
  ChaNGa: x=-3.534189e-09, y=-8.382189e-09, r=9.096790e-09, phi=-1.969810e+00
  delta : dx=1.154636e-12, dy=-1.750873e-09, dr=1.581937e-09, raw_dphi=9.078440e-02, wrapped_dphi=9.078440e-02 (deg=5.202)
Step 38, t=2.387610:
  Quinn : x=-3.736694e-09, y=-6.174401e-09, r=7.217070e-09, phi=-2.115024e+00
  ChaNGa: x=-3.735588e-09, y=-8.583588e-09, r=9.361228e-09, phi=-1.981276e+00
  delta : dx=1.106037e-12, dy=-2.409187e-09, dr=2.144158e-09, raw_dphi=1.337486e-01, wrapped_dphi=1.337486e-01 (deg=7.663)
Step 39, t=2.450442:
  Quinn : x=-3.923293e-09, y=-5.693109e-09, r=6.914023e-09, phi=-2.174190e+00
  ChaNGa: x=-3.922244e-09, y=-8.770244e-09, r=9.607350e-09, phi=-1.991337e+00
  delta : dx=1.048470e-12, dy=-3.077135e-09, dr=2.693327e-09, raw_dphi=1.828522e-01, wrapped_dphi=1.828522e-01 (deg=10.477)
Step 40, t=2.513274:
  Quinn : x=-4.094403e-09, y=-5.189343e-09, r=6.610099e-09, phi=-2.238795e+00
  ChaNGa: x=-4.093421e-09, y=-8.941421e-09, r=9.833876e-09, phi=-2.000121e+00
  delta : dx=9.815507e-13, dy=-3.752079e-09, dr=3.223777e-09, raw_dphi=2.386739e-01, wrapped_dphi=2.386739e-01 (deg=13.675)
Step 41, t=2.576106:
  Quinn : x=-4.249349e-09, y=-4.665089e-09, r=6.310311e-09, phi=-2.309591e+00
  ChaNGa: x=-4.248443e-09, y=-9.096443e-09, r=1.003965e-08, phi=-2.007734e+00
  delta : dx=9.053611e-13, dy=-4.431354e-09, dr=3.729338e-09, raw_dphi=3.018578e-01, wrapped_dphi=3.018578e-01 (deg=17.295)
Step 42, t=2.638938:
  Quinn : x=-4.387519e-09, y=-4.122419e-09, r=6.020354e-09, phi=-2.387336e+00
  ChaNGa: x=-4.386699e-09, y=-9.234698e-09, r=1.022364e-08, phi=-2.014264e+00
  delta : dx=8.202228e-13, dy=-5.112279e-09, dr=4.203284e-09, raw_dphi=3.730724e-01, wrapped_dphi=3.730724e-01 (deg=21.375)
Step 43, t=2.701770:
  Quinn : x=-4.508368e-09, y=-3.563474e-09, r=5.746628e-09, phi=-2.472725e+00
  ChaNGa: x=-4.507642e-09, y=-9.355642e-09, r=1.038493e-08, phi=-2.019786e+00
  delta : dx=7.262246e-13, dy=-5.792168e-09, dr=4.638307e-09, raw_dphi=4.529383e-01, wrapped_dphi=4.529383e-01 (deg=25.951)
Step 44, t=2.764602:
  Quinn : x=-4.611419e-09, y=-2.990461e-09, r=5.496184e-09, phi=-2.566279e+00
  ChaNGa: x=-4.610795e-09, y=-9.458795e-09, r=1.052275e-08, phi=-2.024363e+00
  delta : dx=6.233520e-13, dy=-6.468334e-09, dr=5.026565e-09, raw_dphi=5.419167e-01, wrapped_dphi=5.419167e-01 (deg=31.050)
Step 45, t=2.827433:
  Quinn : x=-4.696264e-09, y=-2.405642e-09, r=5.276553e-09, phi=-2.668196e+00
  ChaNGa: x=-4.695752e-09, y=-9.543752e-09, r=1.063641e-08, phi=-2.028043e+00
  delta : dx=5.122157e-13, dy=-7.138110e-09, dr=5.359861e-09, raw_dphi=6.401539e-01, wrapped_dphi=6.401539e-01 (deg=36.678)
Step 46, t=2.890265:
  Quinn : x=-4.762570e-09, y=-1.811326e-09, r=5.095387e-09, phi=-2.778161e+00
  ChaNGa: x=-4.762177e-09, y=-9.610177e-09, r=1.072538e-08, phi=-2.030865e+00
  delta : dx=3.927674e-13, dy=-7.798851e-09, dr=5.629995e-09, raw_dphi=7.472960e-01, wrapped_dphi=7.472960e-01 (deg=42.817)
Step 47, t=2.953097:
  Quinn : x=-4.810073e-09, y=-1.209859e-09, r=4.959895e-09, phi=-2.895178e+00
  ChaNGa: x=-4.809807e-09, y=-9.657808e-09, r=1.078923e-08, phi=-2.032861e+00
  delta : dx=2.656225e-13, dy=-8.447949e-09, dr=5.829335e-09, raw_dphi=8.623173e-01, wrapped_dphi=8.623173e-01 (deg=49.407)
Step 48, t=3.015929:
  Quinn : x=-4.838587e-09, y=-6.036158e-10, r=4.876092e-09, phi=-3.017483e+00
  ChaNGa: x=-4.838456e-09, y=-9.686456e-09, r=1.082765e-08, phi=-2.034050e+00
  delta : dx=1.310707e-13, dy=-9.082840e-09, dr=5.951561e-09, raw_dphi=9.834337e-01, wrapped_dphi=9.834337e-01 (deg=56.347)
Step 49, t=3.078761:
  Quinn : x=-4.847999e-09, y=5.010366e-12, r=4.848002e-09, phi=3.140559e+00
  ChaNGa: x=-4.848009e-09, y=9.695991e-09, r=1.084045e-08, phi=2.034445e+00
  delta : dx=-9.894737e-15, dy=9.690980e-09, dr=5.992451e-09, raw_dphi=-1.106114e+00, wrapped_dphi=-1.106114e+00 (deg=-63.376)
Step 50, t=3.141593:
  Quinn : x=-4.838272e-09, y=6.136168e-10, r=4.877028e-09, phi=3.015441e+00
  ChaNGa: x=-4.838430e-09, y=-9.686430e-09, r=1.082762e-08, phi=-2.034049e+00
  delta : dx=-1.574442e-13, dy=-1.030005e-08, dr=5.950591e-09, raw_dphi=-5.049489e+00, wrapped_dphi=1.233696e+00 (deg=70.686)
Step 51, t=3.204425:
  Quinn : x=-4.809445e-09, y=1.219801e-09, r=4.961721e-09, phi=2.893204e+00
  ChaNGa: x=-4.809755e-09, y=-9.657755e-09, r=1.078916e-08, phi=-2.032859e+00
  delta : dx=-3.104535e-13, dy=-1.087756e-08, dr=5.827440e-09, raw_dphi=-4.926063e+00, wrapped_dphi=1.357123e+00 (deg=77.757)
Step 52, t=3.267256:
  Quinn : x=-4.761630e-09, y=1.821169e-09, r=5.098017e-09, phi=2.776291e+00
  ChaNGa: x=-4.762099e-09, y=-9.610099e-09, r=1.072528e-08, phi=-2.030862e+00
  delta : dx=-4.686779e-13, dy=-1.143127e-08, dr=5.627261e-09, raw_dphi=-4.807153e+00, wrapped_dphi=1.476032e+00 (deg=84.570)
Step 53, t=3.330088:
  Quinn : x=-4.695017e-09, y=2.415348e-09, r=5.279876e-09, phi=2.666453e+00
  ChaNGa: x=-4.695648e-09, y=-9.543649e-09, r=1.063627e-08, phi=-2.028038e+00
  delta : dx=-6.311790e-13, dy=-1.195900e-08, dr=5.356399e-09, raw_dphi=-4.694491e+00, wrapped_dphi=1.588695e+00 (deg=91.025)
Step 54, t=3.392920:
  Quinn : x=-4.609869e-09, y=2.999991e-09, r=5.500076e-09, phi=2.564672e+00
  ChaNGa: x=-4.610667e-09, y=-9.458666e-09, r=1.052258e-08, phi=-2.024357e+00
  delta : dx=-7.973588e-13, dy=-1.245866e-08, dr=5.022500e-09, raw_dphi=-4.589029e+00, wrapped_dphi=1.694156e+00 (deg=97.068)
Step 55, t=3.455752:
  Quinn : x=-4.506522e-09, y=3.572791e-09, r=5.750963e-09, phi=2.471255e+00
  ChaNGa: x=-4.507489e-09, y=-9.355489e-09, r=1.038473e-08, phi=-2.019779e+00
  delta : dx=-9.662638e-13, dy=-1.292828e-08, dr=4.633767e-09, raw_dphi=-4.491034e+00, wrapped_dphi=1.792151e+00 (deg=102.683)
Step 56, t=3.518584:
  Quinn : x=-4.385384e-09, y=4.131486e-09, r=6.025012e-09, phi=2.385997e+00
  ChaNGa: x=-4.386521e-09, y=-9.234522e-09, r=1.022340e-08, phi=-2.014256e+00
  delta : dx=-1.137046e-12, dy=-1.336601e-08, dr=4.198390e-09, raw_dphi=-4.400252e+00, wrapped_dphi=1.882933e+00 (deg=107.884)
Step 57, t=3.581416:
  Quinn : x=-4.246933e-09, y=4.673870e-09, r=6.315181e-09, phi=2.308372e+00
  ChaNGa: x=-4.248243e-09, y=-9.096243e-09, r=1.003938e-08, phi=-2.007724e+00
  delta : dx=-1.309416e-12, dy=-1.377011e-08, dr=3.724202e-09, raw_dphi=-4.316096e+00, wrapped_dphi=1.967089e+00 (deg=112.706)
Step 58, t=3.644247:
  Quinn : x=-4.091716e-09, y=5.197803e-09, r=6.615081e-09, phi=2.237684e+00
  ChaNGa: x=-4.093198e-09, y=-8.941198e-09, r=9.833580e-09, phi=-2.000110e+00
  delta : dx=-1.481778e-12, dy=-1.413900e-08, dr=3.218499e-09, raw_dphi=-4.237794e+00, wrapped_dphi=2.045391e+00 (deg=117.192)
Step 59, t=3.707079:
  Quinn : x=-3.920346e-09, y=5.701216e-09, r=6.919030e-09, phi=2.173174e+00
  ChaNGa: x=-3.921999e-09, y=-8.770000e-09, r=9.607027e-09, phi=-1.991325e+00
  delta : dx=-1.653671e-12, dy=-1.447122e-08, dr=2.687998e-09, raw_dphi=-4.164499e+00, wrapped_dphi=2.118687e+00 (deg=121.392)
Step 60, t=3.769911:
  Quinn : x=-3.733498e-09, y=6.182121e-09, r=7.222024e-09, phi=2.114092e+00
  ChaNGa: x=-3.735322e-09, y=-8.583322e-09, r=9.360879e-09, phi=-1.981261e+00
  delta : dx=-1.824150e-12, dy=-1.476544e-08, dr=2.138855e-09, raw_dphi=-4.095353e+00, wrapped_dphi=2.187832e+00 (deg=125.354)
Step 61, t=3.832743:
  Quinn : x=-3.531912e-09, y=6.638620e-09, r=7.519686e-09, phi=2.059734e+00
  ChaNGa: x=-3.533904e-09, y=-8.381904e-09, r=9.096417e-09, phi=-1.969793e+00
  delta : dx=-1.992425e-12, dy=-1.502052e-08, dr=1.576731e-09, raw_dphi=-4.029528e+00, wrapped_dphi=2.253658e+00 (deg=129.125)
Step 62, t=3.895575:
  Quinn : x=-3.316382e-09, y=7.068911e-09, r=7.808194e-09, phi=2.009461e+00
  ChaNGa: x=-3.318539e-09, y=-8.166539e-09, r=8.815047e-09, phi=-1.956772e+00
  delta : dx=-2.157025e-12, dy=-1.523545e-08, dr=1.006853e-09, raw_dphi=-3.966233e+00, wrapped_dphi=2.316953e+00 (deg=132.752)
Step 63, t=3.958407:
  Quinn : x=-3.087759e-09, y=7.471295e-09, r=8.084213e-09, phi=1.962701e+00
  ChaNGa: x=-3.090077e-09, y=-7.938077e-09, r=8.518312e-09, phi=-1.942021e+00
  delta : dx=-2.317778e-12, dy=-1.540937e-08, dr=4.340984e-10, raw_dphi=-3.904722e+00, wrapped_dphi=2.378464e+00 (deg=136.276)
Step 64, t=4.021239:
  Quinn : x=-2.846946e-09, y=7.844184e-09, r=8.344838e-09, phi=1.918950e+00
  ChaNGa: x=-2.849419e-09, y=-7.697419e-09, r=8.207890e-09, phi=-1.925333e+00
  delta : dx=-2.473108e-12, dy=-1.554160e-08, dr=-1.369481e-10, raw_dphi=-3.844283e+00, wrapped_dphi=2.438902e+00 (deg=139.739)
Step 65, t=4.084070:
  Quinn : x=-2.594894e-09, y=8.186105e-09, r=8.587537e-09, phi=1.877764e+00
  ChaNGa: x=-2.597517e-09, y=-7.445517e-09, r=7.885608e-09, phi=-1.906464e+00
  delta : dx=-2.622723e-12, dy=-1.563162e-08, dr=-7.019285e-10, raw_dphi=-3.784228e+00, wrapped_dphi=2.498957e+00 (deg=143.180)
Step 66, t=4.146902:
  Quinn : x=-2.332598e-09, y=8.495708e-09, r=8.810112e-09, phi=1.838755e+00
  ChaNGa: x=-2.335363e-09, y=-7.183363e-09, r=7.553451e-09, phi=-1.885125e+00
  delta : dx=-2.765018e-12, dy=-1.567907e-08, dr=-1.256661e-09, raw_dphi=-3.723880e+00, wrapped_dphi=2.559305e+00 (deg=146.637)
Step 67, t=4.209734:
  Quinn : x=-2.061093e-09, y=8.771772e-09, r=9.010665e-09, phi=1.801579e+00
  ChaNGa: x=-2.063993e-09, y=-6.911993e-09, r=7.213578e-09, phi=-1.860978e+00
  delta : dx=-2.899556e-12, dy=-1.568376e-08, dr=-1.797087e-09, raw_dphi=-3.662557e+00, wrapped_dphi=2.620629e+00 (deg=150.151)
Step 68, t=4.272566:
  Quinn : x=-1.781451e-09, y=9.013206e-09, r=9.187570e-09, phi=1.765930e+00
  ChaNGa: x=-1.784477e-09, y=-6.632477e-09, r=6.868341e-09, phi=-1.833624e+00
  delta : dx=-3.025233e-12, dy=-1.564568e-08, dr=-2.319230e-09, raw_dphi=-3.599554e+00, wrapped_dphi=2.683631e+00 (deg=153.761)
Step 69, t=4.335398:
  Quinn : x=-1.494777e-09, y=9.219057e-09, r=9.339453e-09, phi=1.731537e+00
  ChaNGa: x=-1.497918e-09, y=-6.345918e-09, r=6.520310e-09, phi=-1.802598e+00
  delta : dx=-3.141328e-12, dy=-1.556498e-08, dr=-2.819143e-09, raw_dphi=-3.534135e+00, wrapped_dphi=2.749050e+00 (deg=157.509)
Step 70, t=4.398230:
  Quinn : x=-1.202201e-09, y=9.388514e-09, r=9.465172e-09, phi=1.698153e+00
  ChaNGa: x=-1.205448e-09, y=-6.053448e-09, r=6.172304e-09, phi=-1.767359e+00
  delta : dx=-3.247026e-12, dy=-1.544196e-08, dr=-3.292868e-09, raw_dphi=-3.465513e+00, wrapped_dphi=2.817673e+00 (deg=161.441)
Step 71, t=4.461062:
  Quinn : x=-9.048790e-10, y=9.520905e-09, r=9.563809e-09, phi=1.665553e+00
  ChaNGa: x=-9.082204e-10, y=-5.756220e-09, r=5.827430e-09, phi=-1.727287e+00
  delta : dx=-3.341403e-12, dy=-1.527713e-08, dr=-3.736379e-09, raw_dphi=-3.392840e+00, wrapped_dphi=2.890345e+00 (deg=165.605)
Step 72, t=4.523893:
  Quinn : x=-6.039848e-10, y=9.615710e-09, r=9.634660e-09, phi=1.633526e+00
  ChaNGa: x=-6.074080e-10, y=-5.455408e-09, r=5.489118e-09, phi=-1.681680e+00
  delta : dx=-3.423221e-12, dy=-1.507112e-08, dr=-4.145542e-09, raw_dphi=-3.315206e+00, wrapped_dphi=2.967979e+00 (deg=170.053)
Step 73, t=4.586725:
  Quinn : x=-3.007062e-10, y=9.672553e-09, r=9.677227e-09, phi=1.601875e+00
  ChaNGa: x=-3.041984e-10, y=-5.152199e-09, r=5.161171e-09, phi=-1.629770e+00
  delta : dx=-3.492272e-12, dy=-1.482475e-08, dr=-4.516055e-09, raw_dphi=-3.231645e+00, wrapped_dphi=3.051540e+00 (deg=174.840)
Step 74, t=4.649557:
  Quinn : x=3.759630e-12, y=9.691211e-09, r=9.691212e-09, phi=1.570408e+00
  ChaNGa: x=2.134215e-13, y=-4.847786e-09, r=4.847786e-09, phi=-1.570752e+00
  delta : dx=-3.546209e-12, dy=-1.453900e-08, dr=-4.843425e-09, raw_dphi=-3.141161e+00, wrapped_dphi=-3.141161e+00 (deg=-179.975)
Step 75, t=4.712389:
  Quinn : x=3.082106e-10, y=9.671609e-09, r=9.676519e-09, phi=1.538940e+00
  ChaNGa: x=3.048223e-10, y=-4.543178e-09, r=4.553392e-09, phi=-1.503802e+00
  delta : dx=-3.388311e-12, dy=-1.421479e-08, dr=-5.123127e-09, raw_dphi=-3.042742e+00, wrapped_dphi=-3.042742e+00 (deg=-174.336)
Step 76, t=4.775221:
  Quinn : x=6.114448e-10, y=9.613826e-09, r=9.633250e-09, phi=1.507281e+00
  ChaNGa: x=6.080287e-10, y=-4.239971e-09, r=4.283346e-09, phi=-1.428363e+00
  delta : dx=-3.416035e-12, dy=-1.385380e-08, dr=-5.349904e-09, raw_dphi=-2.935645e+00, wrapped_dphi=-2.935645e+00 (deg=-168.200)
Step 77, t=4.838053:
  Quinn : x=9.122651e-10, y=9.518088e-09, r=9.561707e-09, phi=1.475243e+00
  ChaNGa: x=9.088356e-10, y=-3.939164e-09, r=4.042647e-09, phi=-1.344046e+00
  delta : dx=-3.429496e-12, dy=-1.345725e-08, dr=-5.519059e-09, raw_dphi=-2.819289e+00, wrapped_dphi=-2.819289e+00 (deg=-161.533)
Step 78, t=4.900885:
  Quinn : x=1.209484e-09, y=9.384775e-09, r=9.462391e-09, phi=1.442626e+00
  ChaNGa: x=1.206055e-09, y=-3.641945e-09, r=3.836448e-09, phi=-1.251006e+00
  delta : dx=-3.428670e-12, dy=-1.302672e-08, dr=-5.625944e-09, raw_dphi=-2.693631e+00, wrapped_dphi=-2.693631e+00 (deg=-154.334)
Step 79, t=4.963716:
  Quinn : x=1.501928e-09, y=9.214412e-09, r=9.336015e-09, phi=1.409220e+00
  ChaNGa: x=1.498515e-09, y=-3.349485e-09, r=3.669414e-09, phi=-1.150118e+00
  delta : dx=-3.412940e-12, dy=-1.256390e-08, dr=-5.666601e-09, raw_dphi=-2.559337e+00, wrapped_dphi=-2.559337e+00 (deg=-146.639)
Step 80, t=5.026548:
  Quinn : x=1.788442e-09, y=9.007672e-09, r=9.183500e-09, phi=1.374799e+00
  ChaNGa: x=1.785061e-09, y=-3.062939e-09, r=3.545143e-09, phi=-1.043125e+00
  delta : dx=-3.381688e-12, dy=-1.207061e-08, dr=-5.638357e-09, raw_dphi=-2.417923e+00, wrapped_dphi=-2.417923e+00 (deg=-138.537)
Step 81, t=5.089380:
  Quinn : x=2.067897e-09, y=8.765371e-09, r=9.005994e-09, phi=1.339116e+00
  ChaNGa: x=2.064562e-09, y=-2.783438e-09, r=3.465537e-09, phi=-9.326086e-01
  delta : dx=-3.334874e-12, dy=-1.154881e-08, dr=-5.540457e-09, raw_dphi=-2.271724e+00, wrapped_dphi=-2.271724e+00 (deg=-130.160)
Step 82, t=5.152212:
  Quinn : x=2.339187e-09, y=8.488466e-09, r=8.804876e-09, phi=1.301898e+00
  ChaNGa: x=2.335915e-09, y=-2.512085e-09, r=3.430316e-09, phi=-8.217210e-01
  delta : dx=-3.272086e-12, dy=-1.100055e-08, dr=-5.374560e-09, raw_dphi=-2.123619e+00, wrapped_dphi=-2.123619e+00 (deg=-121.674)
Step 83, t=5.215044:
  Quinn : x=2.601243e-09, y=8.178049e-09, r=8.581780e-09, phi=1.262839e+00
  ChaNGa: x=2.598049e-09, y=-2.249951e-09, r=3.436879e-09, phi=-7.137187e-01
  delta : dx=-3.193354e-12, dy=-1.042800e-08, dr=-5.144901e-09, raw_dphi=-1.976558e+00, wrapped_dphi=-1.976558e+00 (deg=-113.248)
Step 84, t=5.277876:
  Quinn : x=2.853029e-09, y=7.835347e-09, r=8.338611e-09, phi=1.221596e+00
  ChaNGa: x=2.849930e-09, y=-1.998070e-09, r=3.480573e-09, phi=-6.114600e-01
  delta : dx=-3.098758e-12, dy=-9.833417e-09, dr=-4.858039e-09, raw_dphi=-1.833056e+00, wrapped_dphi=-1.833056e+00 (deg=-105.026)
Step 85, t=5.340708:
  Quinn : x=3.093552e-09, y=7.461713e-09, r=8.077575e-09, phi=1.177776e+00
  ChaNGa: x=3.090564e-09, y=-1.757436e-09, r=3.555301e-09, phi=-5.170458e-01
  delta : dx=-2.988174e-12, dy=-9.219149e-09, dr=-4.522274e-09, raw_dphi=-1.694822e+00, wrapped_dphi=-1.694822e+00 (deg=-97.106)
Step 86, t=5.403539:
  Quinn : x=3.321862e-09, y=7.058620e-09, r=7.801211e-09, phi=1.130935e+00
  ChaNGa: x=3.319000e-09, y=-1.529000e-09, r=3.654258e-09, phi=-4.317005e-01
  delta : dx=-2.861878e-12, dy=-8.587620e-09, dr=-4.146952e-09, raw_dphi=-1.562636e+00, wrapped_dphi=-1.562636e+00 (deg=-89.532)
Step 87, t=5.466371:
  Quinn : x=3.537058e-09, y=6.627662e-09, r=7.512435e-09, phi=1.080568e+00
  ChaNGa: x=3.534338e-09, y=-1.313662e-09, r=3.770577e-09, phi=-3.558615e-01
  delta : dx=-2.719815e-12, dy=-7.941323e-09, dr=-3.741857e-09, raw_dphi=-1.436430e+00, wrapped_dphi=-1.436430e+00 (deg=-82.301)
Step 88, t=5.529203:
  Quinn : x=3.738290e-09, y=6.170538e-09, r=7.214593e-09, phi=1.026102e+00
  ChaNGa: x=3.735728e-09, y=-1.112272e-09, r=3.897796e-09, phi=-2.893812e-01
  delta : dx=-2.562553e-12, dy=-7.282810e-09, dr=-3.316797e-09, raw_dphi=-1.315483e+00, wrapped_dphi=-1.315483e+00 (deg=-75.372)
Step 89, t=5.592035:
  Quinn : x=3.924765e-09, y=5.689054e-09, r=6.911520e-09, phi=9.668947e-01
  ChaNGa: x=3.922374e-09, y=-9.256258e-10, r=4.030112e-09, phi=-2.317462e-01
  delta : dx=-2.390310e-12, dy=-6.614680e-09, dr=-2.881408e-09, raw_dphi=-1.198641e+00, wrapped_dphi=-1.198641e+00 (deg=-68.677)
Step 90, t=5.654867:
  Quinn : x=4.095744e-09, y=5.185111e-09, r=6.607609e-09, phi=9.022413e-01
  ChaNGa: x=4.093541e-09, y=-7.544594e-10, r=4.162485e-09, phi=-1.822595e-01
  delta : dx=-2.203483e-12, dy=-5.939570e-09, dr=-2.445123e-09, raw_dphi=-1.084501e+00, wrapped_dphi=-1.084501e+00 (deg=-62.137)
Step 91, t=5.717699:
  Quinn : x=4.250555e-09, y=4.660697e-09, r=6.307877e-09, phi=8.313909e-01
  ChaNGa: x=4.248552e-09, y=-5.994483e-10, r=4.290633e-09, phi=-1.401695e-01
  delta : dx=-2.002955e-12, dy=-5.260145e-09, dr=-2.017244e-09, raw_dphi=-9.715604e-01, wrapped_dphi=-9.715604e-01 (deg=-55.666)
Step 92, t=5.780530:
  Quinn : x=4.388584e-09, y=4.117884e-09, r=6.018026e-09, phi=7.535860e-01
  ChaNGa: x=4.386796e-09, y=-4.612043e-10, r=4.410973e-09, phi=-1.047498e-01
  delta : dx=-1.788718e-12, dy=-4.579088e-09, dr=-1.607053e-09, raw_dphi=-8.583358e-01, wrapped_dphi=-8.583358e-01 (deg=-49.179)
Step 93, t=5.843362:
  Quinn : x=4.509289e-09, y=3.558814e-09, r=5.744462e-09, phi=6.681322e-01
  ChaNGa: x=4.507727e-09, y=-3.402729e-10, r=4.520552e-09, phi=-7.534370e-02
  delta : dx=-1.561992e-12, dy=-3.899087e-09, dr=-1.223910e-09, raw_dphi=-7.434759e-01, wrapped_dphi=-7.434759e-01 (deg=-42.598)
Step 94, t=5.906194:
  Quinn : x=4.612191e-09, y=2.985695e-09, r=5.494241e-09, phi=5.745089e-01
  ChaNGa: x=4.610868e-09, y=-2.371316e-10, r=4.616962e-09, phi=-5.138356e-02
  delta : dx=-1.323045e-12, dy=-3.222826e-09, dr=-8.772788e-10, raw_dphi=-6.258925e-01, wrapped_dphi=-6.258925e-01 (deg=-35.861)
Step 95, t=5.969026:
  Quinn : x=4.696886e-09, y=2.400788e-09, r=5.274895e-09, phi=4.725235e-01
  ChaNGa: x=4.695813e-09, y=-1.521872e-10, r=4.698278e-09, phi=-3.239778e-02
  delta : dx=-1.072783e-12, dy=-2.552975e-09, dr=-5.766167e-10, raw_dphi=-5.049213e-01, wrapped_dphi=-5.049213e-01 (deg=-28.930)
Step 96, t=6.031858:
  Quinn : x=4.763037e-09, y=1.806404e-09, r=5.094077e-09, phi=3.624955e-01
  ChaNGa: x=4.762225e-09, y=-8.577501e-11, r=4.762997e-09, phi=-1.800959e-02
  delta : dx=-8.123338e-13, dy=-1.892179e-09, dr=-3.310794e-10, raw_dphi=-3.805051e-01, wrapped_dphi=-3.805051e-01 (deg=-21.801)
Step 97, t=6.094690:
  Quinn : x=4.810385e-09, y=1.204888e-09, r=4.958988e-09, phi=2.454269e-01
  ChaNGa: x=4.809843e-09, y=-3.815718e-11, r=4.809994e-09, phi=-7.932979e-03
  delta : dx=-5.426297e-13, dy=-1.243045e-09, dr=-1.489939e-10, raw_dphi=-2.533599e-01, wrapped_dphi=-2.533599e-01 (deg=-14.516)
Step 98, t=6.157522:
  Quinn : x=4.838743e-09, y=5.986151e-10, r=4.875630e-09, phi=1.230875e-01
  ChaNGa: x=4.838478e-09, y=-9.521600e-12, r=4.838488e-09, phi=-1.967889e-03
  delta : dx=-2.642859e-13, dy=-6.081367e-10, dr=-3.714253e-11, raw_dphi=-1.250554e-01, wrapped_dphi=-1.250554e-01 (deg=-7.165)
Step 99, t=6.220353:
  Quinn : x=4.847997e-09, y=-1.002073e-11, r=4.848008e-09, phi=-2.066981e-03
  ChaNGa: x=4.848018e-09, y=1.871574e-14, r=4.848018e-09, phi=3.860493e-06
  delta : dx=2.108530e-14, dy=1.003945e-11, dr=1.072900e-14, raw_dphi=2.070841e-03, wrapped_dphi=2.070841e-03 (deg=0.119)
ChaNGa trajectory:
  final position (x, y) = (0.00000000, 0.00000000)
Python Quinn integrator:
  nsteps, dt           = 100, 0.062832 (total T = 6.283185, target 6.283185)
  final position (x, y) = (0.00000000, -0.00000000)
Difference at final step (ChaNGa - Quinn):
  dx = 2.109e-14, dy = 1.004e-11, |Δr| = 1.004e-11
Max absolute error over all steps:
  max|dx| = 3.546e-12 kpc, max|dy| = 1.568e-08 kpc
=== Single-particle Hill test (100 steps per orbit, dt = 0.062832) ===
ChaNGa radius : 4.848018e-09 kpc
Quinn radius  : 4.848008e-09 kpc
Δr            : 1.072900e-14 kpc
ChaNGa phase  : 3.860493e-06 rad
Quinn phase   : -2.066981e-03 rad
Δφ            : 2.070841e-03 rad  (0.119 deg)
Saved plot to task2/single_particle_orbit_100steps.png
jinshuozhang@JinshuodeMac-mini Symplectic integrators in the shearing sheet % open task2/single_particle_orbit_100steps.png


