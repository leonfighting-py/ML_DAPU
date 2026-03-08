import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from xgboost import XGBRegressor
import shap

os.chdir('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_1.9')
df = pd.read_csv('Final_Features_Ultimate_all_deletion.csv')

features_mech = [
    'X1_Soft_Content','X2_DA_Content','X3_Hard_Content','X4_R_Ratio','X5_Add_Crosslink',
    'DA_Linker_Type','Network_Topology','Polyol_Class','Iso_Class',
    'Soft_Mw','Constraint_Factor','Hard_Symmetry','Soft_Cryst','Synergy_Feature'
]
features_heal = features_mech + ['healing_temperature','healing_time']

lvl1 = features_mech[:5]
lvl2 = features_mech[5:9]
lvl3 = features_mech[9:] + ['healing_temperature', 'healing_time']

stage_colors = {
    'Stage 1': '#1f77b4',
    'Stage 2': '#2ca02c',
    'Stage 3': '#ff7f0e'
}


def feature_stage(feature_name):
    if feature_name in lvl1:
        return 'Stage 1'
    if feature_name in lvl2:
        return 'Stage 2'
    return 'Stage 3'


def plot_shap_top8_with_rose(target_name, model_name, ranking_df, rose_df, save_path):
    df_bar = ranking_df.sort_values('MeanAbsSHAP', ascending=True).copy()
    bar_colors = [stage_colors[s] for s in df_bar['Stage']]

    fig = plt.figure(figsize=(12, 5), dpi=180)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.7, 1.3], wspace=0.1)

    ax_bar = fig.add_subplot(gs[0, 0])
    bars = ax_bar.barh(
        df_bar['Feature'],
        df_bar['MeanAbsSHAP'],
        color=bar_colors,
        edgecolor='none',
        linewidth=0
    )
    ax_bar.set_xlabel('Mean |SHAP value|')
    ax_bar.set_title(f'{target_name} | {model_name} | SHAP Feature Importance (Top 8)')
    ax_bar.grid(axis='x', linestyle='--', alpha=0.25)
    for spine in ax_bar.spines.values():
        spine.set_linewidth(0.8)
    max_val = float(df_bar['MeanAbsSHAP'].max()) if len(df_bar) else 0.0
    x_pad = max_val * 0.22
    x_lim = max_val + x_pad
    ax_bar.set_xlim(0, x_lim)
    for bar in bars:
        w = bar.get_width()
        x_text = min(w + max_val * 0.015, x_lim - max_val * 0.03)
        ax_bar.text(x_text, bar.get_y() + bar.get_height() / 2, f'{w:.3f}', va='center', fontsize=8)

    legend_handles = [Patch(facecolor=stage_colors[k], edgecolor='none', label=k) for k in ['Stage 1', 'Stage 2', 'Stage 3']]
    ax_bar.legend(handles=legend_handles, loc='lower right', frameon=False)

    ax_rose = fig.add_subplot(gs[0, 1], projection='polar')
    df_rose = rose_df.sort_values('MeanAbsSHAP', ascending=False).copy()
    vals = df_rose['MeanAbsSHAP'].to_numpy(dtype=float)
    pct = vals / vals.sum() * 100
    n = len(df_rose)
    theta = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    width = (2 * np.pi / n) * 0.92
    rose_colors = [stage_colors[s] for s in df_rose['Stage']]

    ax_rose.set_theta_offset(np.pi / 2)
    ax_rose.set_theta_direction(-1)

    primary_n = 4
    r0 = 1.25
    secondary_step = 0.56
    secondary_h = max(pct.max() * 0.07, 0.7)

    for ang in theta:
        ax_rose.plot([ang, ang], [r0, r0 + secondary_step * max(n - primary_n + 1, 1)], color='#e8e8e8', lw=0.7, zorder=0)

    for idx, (ang, share, color) in enumerate(zip(theta, pct, rose_colors)):
        if idx < primary_n:
            bottom = r0
            height = share
        else:
            bottom = r0 + (idx - primary_n + 1) * secondary_step
            height = secondary_h
        ax_rose.bar(
            [ang], [height], width=width,
            bottom=bottom,
            color=color,
            edgecolor='none',
            alpha=0.92,
            zorder=3
        )

    top_radius = max(r0 + pct[:primary_n].max() + 0.7, r0 + secondary_step * (n - primary_n + 2))
    ax_rose.set_ylim(0, top_radius)
    ax_rose.set_xticks([])
    ax_rose.set_yticklabels([])
    ax_rose.grid(alpha=0.18, linewidth=0.6)
    ax_rose.spines['polar'].set_visible(False)

    for idx, (ang, p) in enumerate(zip(theta, pct)):
        if idx < primary_n:
            ax_rose.text(ang, r0 + p + 0.45, f'{p:.1f}%', ha='center', va='center', fontsize=8)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

def get_clean_data(target_col, feature_cols):
    d = df.dropna(subset=[target_col] + feature_cols).copy()
    return d[feature_cols].values, d[target_col].values

def best_model(target_col):
    if target_col == 'tensile_strength':
        return 'SVR', make_pipeline(StandardScaler(), SVR(kernel='rbf', C=100, gamma=0.01, epsilon=0.1))
    if target_col == 'elongation':
        return 'XGBoost', XGBRegressor(
            n_estimators=500, learning_rate=0.02, max_depth=3,
            subsample=0.6, colsample_bytree=0.6, reg_alpha=1.0, reg_lambda=5.0,
            random_state=42, n_jobs=-1
        )
    kernel = 1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-2,1e2), nu=1.5) + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-5,1e1))
    return 'GPR', make_pipeline(StandardScaler(), GaussianProcessRegressor(
        kernel=kernel, alpha=1e-4, normalize_y=False, n_restarts_optimizer=0, random_state=42
    ))

out_dir = 'shap_top8_outputs'
os.makedirs(out_dir, exist_ok=True)

targets = [
    ('Tensile Strength', 'tensile_strength', features_mech),
    ('Elongation', 'elongation', features_mech),
    ('Healing Efficiency', 'healing_eff', features_heal),
]

for target_name, target_col, feats in targets:
    X, y = get_clean_data(target_col, feats)
    X_df = pd.DataFrame(X, columns=feats)
    mname, model = best_model(target_col)
    model.fit(X, y)

    exp = shap.Explainer(model.predict, X_df)(X_df)
    vals = exp.values
    if vals.ndim == 3:
        vals = vals[:, :, 0]

    mean_abs = np.abs(vals).mean(axis=0)
    top_idx = np.argsort(mean_abs)[::-1][:8]
    X_top = X_df.iloc[:, top_idx]
    ranking_df = pd.DataFrame({
        'Feature': X_top.columns,
        'MeanAbsSHAP': mean_abs[top_idx]
    }).sort_values('MeanAbsSHAP', ascending=False)
    ranking_df['Stage'] = ranking_df['Feature'].map(feature_stage)
    ranking_df.to_csv(
        os.path.join(out_dir, f'{target_col}_top8_shap_ranking.csv'), index=False
    )

    all_df = pd.DataFrame({
        'Feature': feats,
        'MeanAbsSHAP': mean_abs
    })
    rose_df = all_df[all_df['Feature'].isin(features_mech)].copy()
    rose_df['Stage'] = rose_df['Feature'].map(feature_stage)

    plot_shap_top8_with_rose(
        target_name=target_name,
        model_name=mname,
        ranking_df=ranking_df,
        rose_df=rose_df,
        save_path=os.path.join(out_dir, f'{target_col}_shap_top8_stage_rose.png')
    )

print('DONE')
print(os.path.abspath(out_dir))
for fn in sorted(os.listdir(out_dir)):
    print(fn)
