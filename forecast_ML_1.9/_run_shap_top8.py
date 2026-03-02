import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    vals_top = vals[:, top_idx]

    pd.DataFrame({
        'Feature': X_top.columns,
        'MeanAbsSHAP': mean_abs[top_idx]
    }).sort_values('MeanAbsSHAP', ascending=False).to_csv(
        os.path.join(out_dir, f'{target_col}_top8_shap_ranking.csv'), index=False
    )

    plt.figure(figsize=(8,5), dpi=140)
    shap.summary_plot(vals_top, X_top, feature_names=X_top.columns.tolist(), plot_type='bar', max_display=8, show=False)
    plt.title(f'{target_name} | {mname} | SHAP Bar (Top 8)')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{target_col}_shap_bar_top8.png'), dpi=300)
    plt.close()

    plt.figure(figsize=(8,5), dpi=140)
    shap.summary_plot(vals_top, X_top, feature_names=X_top.columns.tolist(), max_display=8, show=False)
    plt.title(f'{target_name} | {mname} | SHAP Beeswarm (Top 8)')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{target_col}_shap_beeswarm_top8.png'), dpi=300)
    plt.close()

print('DONE')
print(os.path.abspath(out_dir))
for fn in sorted(os.listdir(out_dir)):
    print(fn)
