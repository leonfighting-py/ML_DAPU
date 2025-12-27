def calculate_recipe(total_mass, x1, x3, r_ratio, 
                     mw_polyol, mw_iso, mw_extender, 
                     f_polyol=2, f_iso=2, f_extender=2):
    """
    根据模型预测的质量分数，计算具体的投料量 (g)
    
    参数:
    total_mass:  总投料量 (g), 例如 30g
    x1:          软段质量分数 (模型输出)
    x3:          硬段质量分数 (模型输出)
    r_ratio:     X4 NCO/OH 比值 (模型输出)
    mw_xxx:      各原料的分子量 (需要查你的试剂瓶)
    f_xxx:       官能度 (默认双官能度)
    """
    
    print(f"--- 投料计算单 (总量 {total_mass}g) ---")
    
    # 1. 计算软段质量 (直接由 X1 决定)
    m_polyol = total_mass * x1
    print(f"1. 多元醇 (Polyol, Mw={mw_polyol}): {m_polyol:.4f} g")
    
    # 2. 计算硬段总质量 (由 X3 决定)
    # 硬段 = 异氰酸酯 + 扩链剂
    m_hard_total = total_mass * x3
    
    # 3. 解方程计算异氰酸酯和扩链剂的分配
    # 方程1 (质量守恒): m_iso + m_ext = m_hard_total
    # 方程2 (官能团守恒): (m_iso/Mw_iso)*f_iso = R * [ (m_polyol/Mw_polyol)*f_polyol + (m_ext/Mw_ext)*f_ext ]
    
    # 这是一个形如 A * m_iso = B + C * m_ext 的方程
    # 推导后求解 m_iso:
    # m_iso = (R * (n_polyol * f_polyol) + R * (m_hard_total / mw_ext) * f_ext) / ( (f_iso / mw_iso) + R * (f_ext / mw_ext) )
    
    n_polyol_groups = (m_polyol / mw_polyol) * f_polyol
    
    numerator = r_ratio * n_polyol_groups + r_ratio * (m_hard_total / mw_extender) * f_extender
    denominator = (f_iso / mw_iso) + r_ratio * (f_extender / mw_extender)
    
    m_iso = numerator / denominator
    m_extender = m_hard_total - m_iso
    
    print(f"2. 异氰酸酯 (Iso, Mw={mw_iso}):   {m_iso:.4f} g")
    print(f"3. 扩链剂 (Extender, Mw={mw_extender}): {m_extender:.4f} g")
    
    # 校验
    print(f"\n[校验] 硬段总质量: {m_iso + m_extender:.4f} g (目标: {m_hard_total:.4f})")
    print(f"[校验] 硬段含量: {(m_iso + m_extender)/total_mass:.2%}")

# === 使用示例 ===
# 假设 Target A 表格里第一行数据是:
# X1=0.45 (PTMG1000), X3=0.55 (MDI), X4(R)=1.02
# 你需要查一下原料的分子量: PTMG=1000, MDI=250.25, BDO(扩链剂)=90.12

calculate_recipe(
    total_mass=30,      # 做 30g 样品
    x1=0.45,            # CSV里的 X1_SoftSeg
    x3=0.55,            # CSV里的 X3_HardSeg
    r_ratio=1.02,       # CSV里的 X4_R_Ratio
    mw_polyol=1000,     # PTMG1000
    mw_iso=250.25,      # MDI
    mw_extender=90.12   # 假设用 BDO 扩链
)