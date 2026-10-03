# 一元弱酸-强碱滴定计算工具

教学实验用滴定曲线计算器：从 JSON 实验配置生成完整滴定曲线，
并核查每个点的物料守恒与电荷平衡。仅使用 Python 标准库（>= 3.12）。

## 模型假设

- 温度固定 25 摄氏度，水的离子积 `Kw = 1e-14`。
- 一元弱酸 `HA <=> H+ + A-`，解离常数 `Ka`；强碱（如 NaOH）完全解离为一价阳离子。
- 忽略活度系数与混合体积收缩，体积直接相加。
- 浓度单位为 mol/L，体积单位为 mL。

每个点同时满足：

```
电荷平衡  [H+] + [Na+] = [OH-] + [A-]
物料守恒  [HA] + [A-] = Ca（稀释后的酸总浓度）
平衡关系  Ka = [H+][A-]/[HA]，Kw = [H+][OH-]
```

对 `h = [H+] > 0` 构造严格单调的电荷失衡函数
`f(h) = h + Cb - Kw/h - Ka*Ca/(Ka+h)`，在 h 的对数域上二分求根
（相对精度 1e-12，上限 200 次迭代，不收敛则报错而非输出伪造结果）。
从初始点、缓冲区、当量点到过量碱区使用完全一致的平衡定义，
不拼接分段近似公式。

## 输入配置

```json
{
  "experiments": [
    {
      "name": "acetic_acid_naoh",
      "acid_volume_ml": 25.0,
      "acid_concentration": 0.1,
      "ka": 1.8e-5,
      "base_concentration": 0.1,
      "base_volumes_ml": [0.0, 12.5, 25.0, 50.0]
    }
  ]
}
```

校验规则：所有数值必须有限；`acid_volume_ml`、`acid_concentration`、`ka`、
`base_concentration` 必须为正；`base_volumes_ml` 允许为零、保留顺序与重复点。
单个实验配置错误会在汇总中报告原因，不影响其他实验。

## 使用方法

```
python -m weak_acid_titration examples/acetic_acid.json --csv out/curve.csv --json out/summary.json
```

- CSV 曲线：每个点的 pH、`[H+]`、`[OH-]`、`[Na+]`、`[A-]`、`[HA]`、
  稀释后浓度、守恒残差与迭代次数。
- JSON 汇总：模型假设、理论当量体积 `Veq = Ca*Va/Cb`、
  每个实验的最大守恒残差与求解迭代上限。

## 验证方法

```
python -m unittest discover -s tests -v
```

测试覆盖：

- 零加碱点与弱酸解析解 `[H+] = (-Ka + sqrt(Ka^2 + 4*Ka*C))/2` 一致；
- 半当量点 pH 接近 pKa（缓冲区）；
- 当量点弱碱性且 `[OH-]` 与 `A-` 水解解析解一致；
- 过量碱区 `[OH-]` 接近过量强碱浓度；
- 全曲线 pH 单调上升，各点电荷/物料守恒残差 < 1e-12 mol/L；
- 批量容错、输入顺序与重复点保留、非法配置拒绝。
