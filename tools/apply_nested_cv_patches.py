#!/usr/bin/env python3
"""
Apply nested LOOCV fixes to:
  - forecast_ML_3.8/final_3.8.ipynb
  - forecast_ML_1.9/alghorithm_multi_1.9_tensile.ipynb
"""
from __future__ import annotations

import json
from pathlib import Path


def to_nb_source(code: str) -> list[str]:
    code = code.rstrip("\n") + "\n"
    return [line + "\n" for line in code.split("\n")[:-1]] + ([code.split("\n")[-1]] if code.split("\n")[-1] else [])


REPO = Path(__file__).resolve().parents[1]

FINAL_CELL6 = r'''from pathlib import Path
import json
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.base import clone
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV, LeaveOneOut, KFold
from sklearn.metrics import r2_score, mean_squared_error

_NB_ROOT = Path.cwd()
if not (_NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv').exists():
    _NB_ROOT = Path.cwd().parent
base_19 = _NB_ROOT / 'forecast_ML_1.9' / 'Final_Features_Ultimate_all_deletion.csv'
base_38 = _NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv'

df_19 = pd.read_csv(base_19)
df_38 = pd.read_csv(base_38)

features_19 = [
    'X1_Soft_Content',
    'X2_DA_Content',
    'X3_Hard_Content',
    'X4_R_Ratio',
    'X5_Add_Crosslink',
    'DA_Linker_Type',
    'Network_Topology',
    'Polyol_Class',
    'Iso_Class',
    'Soft_Mw',
    'Constraint_Factor',
    'Hard_Symmetry',
    'Soft_Cryst',
    'Synergy_Feature',
]

features_38 = [
    'Soft_Segment_Content',
    'Furan_Content',
    'Maleimide_Content',
    'Hard_Segment_Content',
    'NCO_OH_Ratio',
    'Extra_Crosslinker_Content',
    'DA_Linker_Type',
    'Furan_Type',
    'Network_Topology',
    'Polyol_Class',
    'Iso_Class',
    'Soft_Mw',
    'Constraint_Factor',
    'Hard_Symmetry',
    'Soft_Cryst',
    'Synergy_Feature',
]

param_grid = {
    'svr__C': [1, 10, 50, 100, 500, 1000],
    'svr__gamma': ['scale', 0.001, 0.01, 0.1, 0.5, 1.0],
    'svr__epsilon': [0.01, 0.1, 0.5],
}
INNER_RANDOM_STATE = 42


def prepare_xy(df_local, feature_cols, target_col='tensile_strength'):
    data = df_local.dropna(subset=[target_col] + feature_cols).copy()
    X = data[feature_cols]
    y = data[target_col]
    return X, y


def _make_inner_cv(n_train: int):
    n_splits = min(5, n_train)
    if n_splits < 2:
        return None
    return KFold(n_splits=n_splits, shuffle=True, random_state=INNER_RANDOM_STATE)


def _aggregate_params(list_of_dicts: list) -> dict:
    serial = [json.dumps(p, sort_keys=True, default=str) for p in list_of_dicts if p]
    if not serial:
        return {}
    return json.loads(Counter(serial).most_common(1)[0][0])


def run_svr_loocv(label, X, y):
    """外层 LOOCV；内层仅在训练折上做 KFold+GridSearch，避免信息泄露。"""
    print('=' * 100)
    print(f'{label} | SVR 嵌套 LOOCV（内层 KFold 调参）')
    print('=' * 100)
    print(f'样本数: {len(y)} | 特征数: {X.shape[1]}')
    loo = LeaveOneOut()
    y_true = y.to_numpy().ravel()
    y_pred = np.zeros_like(y_true, dtype=float)
    fold_params: list = []
    for tr, te in loo.split(X):
        X_tr, X_te = X.iloc[tr], X.iloc[te]
        y_tr = y.iloc[tr]
        inner = _make_inner_cv(len(tr))
        base = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
        if inner is None:
            m = clone(base)
            m.fit(X_tr, y_tr)
            y_pred[te[0]] = float(m.predict(X_te)[0])
            fold_params.append({})
            continue
        gs = GridSearchCV(
            clone(base),
            param_grid,
            cv=inner,
            scoring='neg_mean_squared_error',
            n_jobs=-1,
            verbose=0,
        )
        gs.fit(X_tr, y_tr)
        y_pred[te[0]] = float(gs.predict(X_te)[0])
        fold_params.append(gs.best_params_)
    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mode_params = _aggregate_params(fold_params)
    print('\n[各折参数众数 | 仅作报告/全数据再训练参考]')
    print(mode_params)
    print(f'R² Score: {r2:.4f}')
    print(f'RMSE: {rmse:.4f}')
    return {
        'Dataset': label,
        'Samples': len(y),
        'Features': X.shape[1],
        'R2': r2,
        'RMSE': rmse,
        'Best_Params': mode_params,
        'y_true': y_true,
        'y_pred': y_pred,
    }

X_19, y_19 = prepare_xy(df_19, features_19)
X_38, y_38 = prepare_xy(df_38, features_38)
result_19 = run_svr_loocv('1.9 原始 tensile 特征集', X_19, y_19)
result_38 = run_svr_loocv('3.8 当前 tensile 特征集', X_38, y_38)
comparison_df = pd.DataFrame([
    {k: v for k, v in result_19.items() if k not in {'Best_Params', 'y_true', 'y_pred'}},
    {k: v for k, v in result_38.items() if k not in {'Best_Params', 'y_true', 'y_pred'}},
])
r2_delta = result_38['R2'] - result_19['R2']
rmse_delta = result_38['RMSE'] - result_19['RMSE']
print('\n' + '=' * 100)
print('SVR 对比汇总')
print('=' * 100)
display(comparison_df)
print(f"R² 差异 (3.8 - 1.9): {r2_delta:+.4f}")
print(f"RMSE 差异 (3.8 - 1.9): {rmse_delta:+.4f}")
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), dpi=120)
axes[0].bar(comparison_df['Dataset'], comparison_df['R2'], color=['#9ecae1', '#F28E2B'], edgecolor='black')
axes[0].set_title('SVR LOOCV R² 对比', fontweight='bold')
axes[0].set_ylabel('R²')
axes[0].tick_params(axis='x', rotation=15)
for idx, val in enumerate(comparison_df['R2']):
    axes[0].text(idx, val + 0.01, f'{val:.3f}', ha='center', va='bottom')
axes[1].bar(comparison_df['Dataset'], comparison_df['RMSE'], color=['#9ecae1', '#F28E2B'], edgecolor='black')
axes[1].set_title('SVR LOOCV RMSE 对比', fontweight='bold')
axes[1].set_ylabel('RMSE')
axes[1].tick_params(axis='x', rotation=15)
for idx, val in enumerate(comparison_df['RMSE']):
    axes[1].text(idx, val + max(comparison_df['RMSE']) * 0.02, f'{val:.2f}', ha='center', va='bottom')
plt.tight_layout()
plt.show()
'''

NESTED_BLOCK_10 = r'''
from sklearn.base import clone
from sklearn.model_selection import ParameterGrid, RandomizedSearchCV


def _inner_kfold_10(n_train: int, random_state: int = 42):
    ns = min(5, n_train)
    if ns < 2:
        return None
    return KFold(n_splits=ns, shuffle=True, random_state=random_state)


def _maj_params(rows: list) -> dict:
    import json as _j
    from collections import Counter as _C
    keys = [_j.dumps(p, sort_keys=True, default=str) for p in rows if p]
    if not keys:
        return {}
    return _j.loads(_C(keys).most_common(1)[0][0])


def nested_loocv_search_eval(model_name, base_estimator, param_space, X, y):
    """外层 LeaveOneOut；调参仅在训练折内完成（内层 KFold）。"""
    loo = LeaveOneOut()
    y = np.asarray(y).ravel()
    y_pred = np.zeros(len(y))
    inner_rmse = []
    fold_params: list = []
    route = 'NestedLOOCV'
    if isinstance(X, pd.DataFrame):
        X_df = X
    else:
        X_df = pd.DataFrame(X)
    y_ser = pd.Series(y)
    for tr, te in loo.split(X_df):
        X_tr, X_te = X_df.iloc[tr], X_df.iloc[te]
        y_tr = y_ser.iloc[tr]
        inner = _inner_kfold_10(len(tr))
        est = clone(base_estimator)
        if not param_space:
            est.fit(X_tr, y_tr)
            y_pred[te[0]] = float(est.predict(X_te)[0])
            fold_params.append({})
            continue
        if inner is None:
            est.fit(X_tr, y_tr)
            y_pred[te[0]] = float(est.predict(X_te)[0])
            fold_params.append({})
            continue
        n_comb = len(list(ParameterGrid(param_space)))
        if n_comb <= GRID_THRESHOLD:
            search = GridSearchCV(est, param_space, cv=inner, scoring='neg_mean_squared_error', n_jobs=-1, verbose=0)
            route = f'NestedGrid({n_comb})'
        else:
            ni = min(RANDOM_ITER, n_comb)
            search = RandomizedSearchCV(
                est, param_distributions=param_space, n_iter=ni, cv=inner,
                scoring='neg_mean_squared_error', random_state=42, n_jobs=-1, verbose=0,
            )
            route = f'NestedRandom({ni}/{n_comb})'
        search.fit(X_tr, y_tr)
        y_pred[te[0]] = float(search.predict(X_te)[0])
        inner_rmse.append(float(np.sqrt(-search.best_score_)))
        fold_params.append(search.best_params_)
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    tune_rmse = float(np.mean(inner_rmse)) if inner_rmse else np.nan
    maj = _maj_params(fold_params)
    final_m = clone(base_estimator)
    if maj:
        final_m.set_params(**maj)
    final_m.fit(X_df, y_ser)
    return final_m, route, maj, tune_rmse, r2, rmse, y_pred

'''


def patch_final(path: Path) -> None:
    nb = json.loads(path.read_text(encoding="utf-8"))

    nb["cells"][6]["source"] = to_nb_source(FINAL_CELL6)

    s10 = "".join(nb["cells"][10]["source"])
    s10 = s10.replace(
        "file_path_38 = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_3.8/Final_Features_3.8.csv')",
        "_NB_ROOT = Path.cwd()\nif not (_NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv').exists():\n    _NB_ROOT = Path.cwd().parent\nfile_path_38 = _NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv'",
    )
    if "nested_loocv_search_eval" not in s10:
        s10 = s10.replace("RANDOM_ITER = 50\n\nbase_models = {", "RANDOM_ITER = 50\n" + NESTED_BLOCK_10 + "\nbase_models = {")
    s10 = s10.replace(
        "def evaluate_model_loocv(model, X, y):\n    y_pred = cross_val_predict(model, X, y, cv=LeaveOneOut(), n_jobs=-1)\n",
        "def evaluate_model_loocv(model, X, y):\n    \"\"\"无调参模型使用；有调参模型请用 nested_loocv_search_eval。\"\"\"\n    y_pred = cross_val_predict(model, X, y, cv=LeaveOneOut(), n_jobs=-1)\n",
    )
    s10 = s10.replace(
        "    tuned_model, route_used, best_params, tune_rmse = tune_model_with_dynamic_route(\n        model_name, base_models[model_name], param_spaces[model_name], X_model, y_model\n    )\n    r2, rmse, _ = evaluate_model_loocv(tuned_model, X_model, y_model)\n",
        "    tuned_model, route_used, best_params, tune_rmse, r2, rmse, _ = nested_loocv_search_eval(\n        model_name, base_models[model_name], param_spaces[model_name], X_model, y_model\n    )\n",
    )
    nb["cells"][10]["source"] = to_nb_source(s10)

    s12 = "".join(nb["cells"][12]["source"])
    s12 = s12.replace(
        "from sklearn.model_selection import LeaveOneOut, cross_val_predict, GridSearchCV\n",
        "from sklearn.base import clone\nfrom sklearn.model_selection import LeaveOneOut, cross_val_predict, GridSearchCV, KFold\n",
    )
    s12 = s12.replace(
        "file_path_strict = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_3.8/Final_Features_3.8.csv')",
        "_NB_ROOT = Path.cwd()\nif not (_NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv').exists():\n    _NB_ROOT = Path.cwd().parent\nfile_path_strict = _NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv'",
    )
    OLD_STRICT_FN = """def strict_tune_model(model_name, model, param_space, X, y):
    if not param_space:
        model.fit(X, y)
        train_pred = model.predict(X)
        train_rmse = np.sqrt(mean_squared_error(y, train_pred))
        return model, 'Default(NoTune)', {}, train_rmse

    total_combinations = 1
    for values in param_space.values():
        total_combinations *= len(values)

    print(f'  -> 组合数: {total_combinations}')
    search = GridSearchCV(
        model,
        param_space,
        cv=strict_cv,
        scoring='neg_mean_squared_error',
        n_jobs=-1,
    )
    search.fit(X, y)
    best_model = search.best_estimator_
    best_params = search.best_params_
    best_rmse_cv = np.sqrt(-search.best_score_)
    return best_model, f'StrictGridSearch({total_combinations}组合, LOOCV)', best_params, best_rmse_cv

"""
    NEW_STRICT_FN = """def strict_tune_model(model_name, model, param_space, X, y):
    import json as _j
    from collections import Counter as _C
    if not param_space:
        r2, rmse, _ = evaluate_loocv(model, X, y)
        model.fit(X, y)
        train_rmse = np.sqrt(mean_squared_error(y, model.predict(X)))
        return model, 'Default(NoTune)', {}, train_rmse, r2, rmse
    total_combinations = 1
    for values in param_space.values():
        total_combinations *= len(values)
    print(f'  -> 组合数: {total_combinations} (内层 KFold 调参)')
    yv = np.asarray(y).ravel()
    y_pred = np.zeros(len(yv))
    inner_scores = []
    fold_params: list = []
    Xdf = X if isinstance(X, pd.DataFrame) else pd.DataFrame(X)
    yser = pd.Series(yv)
    for tr, te in strict_cv.split(Xdf):
        X_tr, X_te = Xdf.iloc[tr], Xdf.iloc[te]
        y_tr = yser.iloc[tr]
        nn = len(tr)
        ns = min(5, nn)
        inner = KFold(n_splits=ns, shuffle=True, random_state=42) if ns >= 2 else None
        est = clone(model)
        if inner is None:
            est.fit(X_tr, y_tr)
            y_pred[te[0]] = float(est.predict(X_te)[0])
            fold_params.append({})
            continue
        search = GridSearchCV(
            est, param_space, cv=inner, scoring='neg_mean_squared_error', n_jobs=-1, verbose=0,
        )
        search.fit(X_tr, y_tr)
        y_pred[te[0]] = float(search.predict(X_te)[0])
        inner_scores.append(float(np.sqrt(-search.best_score_)))
        fold_params.append(search.best_params_)
    r2 = r2_score(yv, y_pred)
    rmse = np.sqrt(mean_squared_error(yv, y_pred))
    tune_rmse = float(np.mean(inner_scores)) if inner_scores else np.nan
    keys = [_j.dumps(p, sort_keys=True, default=str) for p in fold_params if p]
    maj = _j.loads(_C(keys).most_common(1)[0][0]) if keys else {}
    final_m = clone(model)
    if maj:
        final_m.set_params(**maj)
    final_m.fit(Xdf, yser)
    route = f'StrictNestedGrid({total_combinations}组合)'
    return final_m, route, maj, tune_rmse, r2, rmse

"""
    if OLD_STRICT_FN not in s12:
        raise RuntimeError("final_3.8 cell12: expected strict_tune_model block missing")
    s12 = s12.replace(OLD_STRICT_FN, NEW_STRICT_FN)
    s12 = s12.replace(
        "    tuned_model, route_used, best_params, tune_rmse = strict_tune_model(\n        model_name, model, strict_param_spaces[model_name], X_strict, y_strict\n    )\n    r2, rmse, _ = evaluate_loocv(tuned_model, X_strict, y_strict)\n",
        "    tuned_model, route_used, best_params, tune_rmse, r2, rmse = strict_tune_model(\n        model_name, model, strict_param_spaces[model_name], X_strict, y_strict\n    )\n",
    )
    nb["cells"][12]["source"] = to_nb_source(s12)

    s14 = "".join(nb["cells"][14]["source"])
    s14 = s14.replace(
        "ablation_file = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_3.8/Final_Features_3.8.csv')",
        "_NB_ROOT = Path.cwd()\nif not (_NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv').exists():\n    _NB_ROOT = Path.cwd().parent\nablation_file = _NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv'",
    )
    s14 = s14.replace(
        "from sklearn.model_selection import LeaveOneOut, GridSearchCV, cross_val_predict\n",
        "from sklearn.base import clone\nfrom sklearn.model_selection import LeaveOneOut, GridSearchCV, cross_val_predict, KFold\n",
    )
    OLD_SVR = """def tune_and_eval_svr(X, y, param_grid):
    model = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
    search = GridSearchCV(
        model,
        param_grid,
        cv=loo,
        scoring='neg_mean_squared_error',
        n_jobs=-1,
    )
    search.fit(X, y)
    best_model = search.best_estimator_
    y_pred = cross_val_predict(best_model, X, y, cv=loo, n_jobs=-1)
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    return {
        'best_model': best_model,
        'best_params': search.best_params_,
        'r2': r2,
        'rmse': rmse,
    }

"""
    NEW_SVR = """def tune_and_eval_svr(X, y, param_grid):
    yv = np.asarray(y).ravel()
    y_pred = np.zeros(len(yv))
    fold_best = []
    Xdf = X if isinstance(X, pd.DataFrame) else pd.DataFrame(X)
    yser = pd.Series(yv)
    for tr, te in loo.split(Xdf):
        X_tr, X_te = Xdf.iloc[tr], Xdf.iloc[te]
        y_tr = yser.iloc[tr]
        nn = len(tr)
        ns = min(5, nn)
        inner = KFold(n_splits=ns, shuffle=True, random_state=42) if ns >= 2 else None
        base = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
        if inner is None:
            m = clone(base)
            m.fit(X_tr, y_tr)
            y_pred[te[0]] = float(m.predict(X_te)[0])
            fold_best.append({})
            continue
        search = GridSearchCV(
            clone(base), param_grid, cv=inner, scoring='neg_mean_squared_error', n_jobs=-1, verbose=0,
        )
        search.fit(X_tr, y_tr)
        y_pred[te[0]] = float(search.predict(X_te)[0])
        fold_best.append(search.best_params_)
    r2 = r2_score(yv, y_pred)
    rmse = np.sqrt(mean_squared_error(yv, y_pred))
    import json as _j
    from collections import Counter as _C
    keys = [_j.dumps(p, sort_keys=True, default=str) for p in fold_best if p]
    maj = _j.loads(_C(keys).most_common(1)[0][0]) if keys else {}
    best_model = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
    if maj:
        best_model.set_params(**maj)
    best_model.fit(Xdf, yser)
    return {'best_model': best_model, 'best_params': maj, 'r2': r2, 'rmse': rmse}

"""
    if OLD_SVR not in s14:
        raise RuntimeError("final_3.8 cell14: tune_and_eval_svr not found")
    s14 = s14.replace(OLD_SVR, NEW_SVR)
    nb["cells"][14]["source"] = to_nb_source(s14)

    s16 = "".join(nb["cells"][16]["source"])
    s16 = s16.replace(
        "compare_file = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_3.8/Final_Features_3.8.csv')",
        "_NB_ROOT = Path.cwd()\nif not (_NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv').exists():\n    _NB_ROOT = Path.cwd().parent\ncompare_file = _NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv'",
    )
    s16 = s16.replace(
        "old_file = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_1.9/Final_Features_Ultimate_all_deletion.csv')",
        "old_file = _NB_ROOT / 'forecast_ML_1.9' / 'Final_Features_Ultimate_all_deletion.csv'",
    )
    s16 = s16.replace(
        "from sklearn.model_selection import LeaveOneOut, GridSearchCV, cross_val_predict\n",
        "from sklearn.base import clone\nfrom sklearn.model_selection import LeaveOneOut, GridSearchCV, cross_val_predict, KFold\n",
    )
    OLD_ENC = """def run_encoding_svr(X, y):
    model = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
    search = GridSearchCV(
        model,
        encoding_grid,
        cv=encoding_loo,
        scoring='neg_mean_squared_error',
        n_jobs=-1,
    )
    search.fit(X, y)
    best_model = search.best_estimator_
    y_pred = cross_val_predict(best_model, X, y, cv=encoding_loo, n_jobs=-1)
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    return {
        'R2': r2,
        'RMSE': rmse,
        'Best_Params': search.best_params_,
        'Feature_Count': X.shape[1],
    }

"""
    NEW_ENC = """def run_encoding_svr(X, y):
    yv = np.asarray(y).ravel()
    y_pred = np.zeros(len(yv))
    fold_best = []
    Xdf = X.copy() if isinstance(X, pd.DataFrame) else pd.DataFrame(X)
    yser = pd.Series(yv)
    for tr, te in encoding_loo.split(Xdf):
        X_tr, X_te = Xdf.iloc[tr], Xdf.iloc[te]
        y_tr = yser.iloc[tr]
        nn = len(tr)
        ns = min(5, nn)
        inner = KFold(n_splits=ns, shuffle=True, random_state=42) if ns >= 2 else None
        base = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
        if inner is None:
            m = clone(base)
            m.fit(X_tr, y_tr)
            y_pred[te[0]] = float(m.predict(X_te)[0])
            fold_best.append({})
            continue
        search = GridSearchCV(
            clone(base), encoding_grid, cv=inner, scoring='neg_mean_squared_error', n_jobs=-1, verbose=0,
        )
        search.fit(X_tr, y_tr)
        y_pred[te[0]] = float(search.predict(X_te)[0])
        fold_best.append(search.best_params_)
    r2 = r2_score(yv, y_pred)
    rmse = np.sqrt(mean_squared_error(yv, y_pred))
    import json as _j
    from collections import Counter as _C
    keys = [_j.dumps(p, sort_keys=True, default=str) for p in fold_best if p]
    maj = _j.loads(_C(keys).most_common(1)[0][0]) if keys else {}
    return {'R2': r2, 'RMSE': rmse, 'Best_Params': maj, 'Feature_Count': X.shape[1]}

"""
    if OLD_ENC not in s16:
        raise RuntimeError("final_3.8 cell16: run_encoding_svr not found")
    s16 = s16.replace(OLD_ENC, NEW_ENC)
    nb["cells"][16]["source"] = to_nb_source(s16)

    s18 = "".join(nb["cells"][18]["source"])
    s18 = s18.replace(
        """features_mech_38_replica = [
    'Soft_Segment_Content',
    'Furan_Content',
    'Maleimide_Content',
    'Hard_Segment_Content',
    'NCO_OH_Ratio',
    'Extra_Crosslinker_Content',
    'Furan_Type',
    'Network_Topology',
    'Add_Crosslinker_Flag',
    'Polyol_Class',
    'Iso_Class',
    'Soft_Mw',
    'Constraint_Factor',
    'Hard_Symmetry',
    'Soft_Cryst',
    'Synergy_Feature',
    'Extender_Class',
]
""",
        """# 与主流程 3.8 特征口径一致（含 DA_Linker_Type；不含高共线/已弃用列）
features_mech_38_replica = [
    'Soft_Segment_Content',
    'Furan_Content',
    'Maleimide_Content',
    'Hard_Segment_Content',
    'NCO_OH_Ratio',
    'Extra_Crosslinker_Content',
    'DA_Linker_Type',
    'Furan_Type',
    'Network_Topology',
    'Polyol_Class',
    'Iso_Class',
    'Soft_Mw',
    'Constraint_Factor',
    'Hard_Symmetry',
    'Soft_Cryst',
    'Synergy_Feature',
]
""",
    )
    s18 = s18.replace(
        "replica_file = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_3.8/Final_Features_3.8.csv')",
        "_NB_ROOT = Path.cwd()\nif not (_NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv').exists():\n    _NB_ROOT = Path.cwd().parent\nreplica_file = _NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv'",
    )
    nb["cells"][18]["source"] = to_nb_source(s18)

    s20 = "".join(nb["cells"][20]["source"])
    s20 = s20.replace(
        "new_file = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_3.8/Final_Features_3.8.csv')",
        "_NB_ROOT = Path.cwd()\nif not (_NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv').exists():\n    _NB_ROOT = Path.cwd().parent\nnew_file = _NB_ROOT / 'forecast_ML_3.8' / 'Final_Features_3.8.csv'",
    )
    s20 = s20.replace(
        "old_file = Path('/Users/yljj0528/Documents/GitHub/ML_DAPU/forecast_ML_1.9/Final_Features_Ultimate_all_deletion.csv')",
        "old_file = _NB_ROOT / 'forecast_ML_1.9' / 'Final_Features_Ultimate_all_deletion.csv'",
    )
    s20 = s20.replace(
        "from sklearn.model_selection import LeaveOneOut, GridSearchCV, cross_val_predict\n",
        "from sklearn.base import clone\nfrom sklearn.model_selection import LeaveOneOut, GridSearchCV, cross_val_predict, KFold\n",
    )
    OLD_MULTI = """def strict_tune_and_eval_multi(model_name, model, param_space, X, y):
    if not param_space:
        model.fit(X, y)
        y_pred = cross_val_predict(model, X, y, cv=strict_cv_multi, n_jobs=-1)
        r2 = r2_score(y, y_pred)
        rmse = np.sqrt(mean_squared_error(y, y_pred))
        return {
            'best_model': model,
            'route': 'Default(NoTune)',
            'best_params': {},
            'r2': r2,
            'rmse': rmse,
        }

    total_combinations = 1
    for values in param_space.values():
        total_combinations *= len(values)

    search = GridSearchCV(
        model,
        param_space,
        cv=strict_cv_multi,
        scoring='neg_mean_squared_error',
        n_jobs=-1,
    )
    search.fit(X, y)
    best_model = search.best_estimator_
    y_pred = cross_val_predict(best_model, X, y, cv=strict_cv_multi, n_jobs=-1)
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    return {
        'best_model': best_model,
        'route': f'StrictGridSearch({total_combinations}组合, LOOCV)',
        'best_params': search.best_params_,
        'r2': r2,
        'rmse': rmse,
    }

"""
    NEW_MULTI = """def strict_tune_and_eval_multi(model_name, model, param_space, X, y):
    if not param_space:
        model.fit(X, y)
        y_pred = cross_val_predict(model, X, y, cv=strict_cv_multi, n_jobs=-1)
        r2 = r2_score(y, y_pred)
        rmse = np.sqrt(mean_squared_error(y, y_pred))
        return {
            'best_model': model,
            'route': 'Default(NoTune)',
            'best_params': {},
            'r2': r2,
            'rmse': rmse,
        }
    total_combinations = 1
    for values in param_space.values():
        total_combinations *= len(values)
    yv = np.asarray(y).ravel()
    y_pred = np.zeros(len(yv))
    inner_scores = []
    fold_params: list = []
    Xdf = X if isinstance(X, pd.DataFrame) else pd.DataFrame(X)
    yser = pd.Series(yv)
    for tr, te in strict_cv_multi.split(Xdf):
        X_tr, X_te = Xdf.iloc[tr], Xdf.iloc[te]
        y_tr = yser.iloc[tr]
        nn = len(tr)
        ns = min(5, nn)
        inner = KFold(n_splits=ns, shuffle=True, random_state=42) if ns >= 2 else None
        est = clone(model)
        if inner is None:
            est.fit(X_tr, y_tr)
            y_pred[te[0]] = float(est.predict(X_te)[0])
            fold_params.append({})
            continue
        search = GridSearchCV(
            est, param_space, cv=inner, scoring='neg_mean_squared_error', n_jobs=-1, verbose=0,
        )
        search.fit(X_tr, y_tr)
        y_pred[te[0]] = float(search.predict(X_te)[0])
        inner_scores.append(float(np.sqrt(-search.best_score_)))
        fold_params.append(search.best_params_)
    r2 = r2_score(yv, y_pred)
    rmse = np.sqrt(mean_squared_error(yv, y_pred))
    import json as _j
    from collections import Counter as _C
    keys = [_j.dumps(p, sort_keys=True, default=str) for p in fold_params if p]
    maj = _j.loads(_C(keys).most_common(1)[0][0]) if keys else {}
    best_model = clone(model)
    if maj:
        best_model.set_params(**maj)
    best_model.fit(Xdf, yser)
    return {
        'best_model': best_model,
        'route': f'StrictNestedGrid({total_combinations}组合)',
        'best_params': maj,
        'r2': r2,
        'rmse': rmse,
    }

"""
    if OLD_MULTI not in s20:
        raise RuntimeError("final_3.8 cell20: strict_tune_and_eval_multi not found")
    s20 = s20.replace(OLD_MULTI, NEW_MULTI)
    nb["cells"][20]["source"] = to_nb_source(s20)

    path.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print("Patched", path)


TENSILE_CELL8 = r'''# SVR — 嵌套 LOOCV（外层留一，内层 KFold 调参；避免全数据调参后再 LOOCV）
from sklearn.base import clone
from sklearn.model_selection import LeaveOneOut, KFold, GridSearchCV
import json
from collections import Counter

target_col = 'tensile_strength'
features = features_mech
X, y = get_clean_data(target_col, features, df)

param_grid = {
    'svr__C': [1, 10, 50, 100, 500, 1000],
    'svr__gamma': ['scale', 0.001, 0.01, 0.1, 0.5, 1.0],
    'svr__epsilon': [0.01, 0.1, 0.5],
}

def _svr_inner_cv(n_tr: int):
    ns = min(5, n_tr)
    if ns < 2:
        return None
    return KFold(n_splits=ns, shuffle=True, random_state=42)

loo = LeaveOneOut()
y_arr = np.asarray(y).ravel()
X_arr = np.asarray(X)
y_pred_svr = np.zeros(len(y_arr))
fold_params: list = []

for tr, te in loo.split(X_arr):
    X_tr, X_te = X_arr[tr], X_arr[te]
    y_tr = y_arr[tr]
    inner = _svr_inner_cv(len(tr))
    pipe = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
    if inner is None:
        m = clone(pipe)
        m.fit(X_tr, y_tr)
        y_pred_svr[te[0]] = float(m.predict(X_te)[0])
        fold_params.append({})
        continue
    gs = GridSearchCV(
        clone(pipe), param_grid, cv=inner,
        scoring='neg_mean_squared_error', n_jobs=-1, verbose=0,
    )
    gs.fit(X_tr, y_tr)
    y_pred_svr[te[0]] = float(gs.predict(X_te)[0])
    fold_params.append(gs.best_params_)

r2 = r2_score(y_arr, y_pred_svr)
rmse = np.sqrt(mean_squared_error(y_arr, y_pred_svr))

_keys = [json.dumps(p, sort_keys=True, default=str) for p in fold_params if p]
best_params_svr_nested = json.loads(Counter(_keys).most_common(1)[0][0]) if _keys else {}

best_svr = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
if best_params_svr_nested:
    best_svr.set_params(**best_params_svr_nested)
best_svr.fit(X_arr, y_arr)

print("\n[各折参数众数]", best_params_svr_nested)
print(f"\n=== SVR (嵌套 LOOCV) ===\nR2 Score: {r2:.4f}\nRMSE: {rmse:.4f}")

plt.figure(figsize=(6, 5))
min_val = min(y_arr.min(), y_pred_svr.min())
max_val = max(y_arr.max(), y_pred_svr.max())
x_line = np.linspace(min_val, max_val, 200)
plt.fill_between(x_line, x_line - rmse, x_line + rmse, color='gray', alpha=0.2)
plt.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label='Perfect Fit')
plt.plot(x_line, x_line + rmse, 'k--', linewidth=1)
plt.plot(x_line, x_line - rmse, 'k--', linewidth=1)
plt.scatter(y_arr, y_pred_svr, alpha=0.7, c='lightyellow', edgecolors='k', label='Prediction', s=80)
plt.xlabel(f'Measured {target_col}')
plt.ylabel(f'Predicted {target_col}')
plt.title(f'SVR Prediction (nested LOOCV)\nR2={r2:.2f}, RMSE={rmse:.2f}')
plt.legend()
plt.tight_layout()
plt.show()
'''


TENSILE_CELL9 = r'''# SVR + SHAP / Permutation：在独立留出集上算重要性（缓解训练集乐观偏差）
from sklearn.model_selection import train_test_split

X_arr = np.asarray(X)
y_arr = np.asarray(y).ravel()
X_pi_tr, X_pi_val, y_pi_tr, y_pi_val = train_test_split(
    X_arr, y_arr, test_size=0.2, random_state=42, shuffle=True
)

interp_svr = make_pipeline(StandardScaler(), SVR(kernel='rbf'))
_params = globals().get('best_params_svr_nested', {})
if _params:
    interp_svr.set_params(**_params)
interp_svr.fit(X_pi_tr, y_pi_tr)

print("\n正在计算排列重要性 (Permutation Importance, 20% 留出集)...")
r = permutation_importance(
    interp_svr, X_pi_val, y_pi_val, n_repeats=30, random_state=42, scoring='r2'
)
importances = r.importances_mean
indices = np.argsort(importances)[::-1]

print("\n[特征重要性排行 (基于 R2 损失)]")
top_features_perm = []
for i in range(len(features)):
    idx = indices[i]
    if importances[idx] > 0:
        print(f"{i+1}. {features[idx]}: {importances[idx]:.4f}")
        top_features_perm.append((features[idx], importances[idx]))

print("\n正在计算 SVR 的 SHAP 特征重要性 ...")
X_val_sub = X_pi_val[: min(50, len(X_pi_val))]
X_bg = shap.kmeans(X_pi_tr, k=min(30, len(X_pi_tr))).data
explainer = shap.KernelExplainer(interp_svr.predict, X_bg)
shap_values = explainer.shap_values(X_val_sub, nsamples=min(500, 2000))
shap_importances = np.abs(np.asarray(shap_values)).mean(axis=0)
shap_indices = np.argsort(shap_importances)[::-1]

print("\n[SHAP 特征重要性排行 (基于平均绝对贡献, 留出子集)]")
top_features_shap = []
for i in range(len(features)):
    idx = shap_indices[i]
    print(f"{i+1}. {features[idx]}: {shap_importances[idx]:.4f}")
    top_features_shap.append((features[idx], shap_importances[idx]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
ax1.bar(range(len(top_features_perm)), [val for _, val in top_features_perm], color='#4292c6', align="center")
ax1.set_title("SVR Feature Importance (Permutation, held-out)", fontsize=14)
ax1.set_xticks(range(len(top_features_perm)))
ax1.set_xticklabels([name for name, _ in top_features_perm], rotation=45, ha='right')
ax1.set_ylabel("Decrease in R2 Score")

ax2.bar(range(len(top_features_shap)), [val for _, val in top_features_shap], color='lightgreen', align="center")
ax2.set_title("SVR Feature Importance (SHAP, held-out subset)", fontsize=14)
ax2.set_xticks(range(len(top_features_shap)))
ax2.set_xticklabels([name for name, _ in top_features_shap], rotation=45, ha='right')
ax2.set_ylabel("Mean |SHAP Value|")

plt.tight_layout()
plt.show()

shap.summary_plot(shap_values, X_val_sub, feature_names=features, show=False)
plt.title(f"SVR SHAP Summary Plot - Target: {target_col}")
plt.tight_layout()
plt.show()
'''


def patch_tensile(path: Path) -> None:
    nb = json.loads(path.read_text(encoding="utf-8"))

    s3 = "".join(nb["cells"][3]["source"])
    NESTED_FN = r'''

# 6. 嵌套 LOOCV + 内层 KFold 调参（多模型单元调用）
from sklearn.base import clone
from sklearn.model_selection import ParameterGrid, RandomizedSearchCV


def nested_loocv_tune_eval(model_name, base_estimator, param_space, X, y, grid_threshold=60, random_iter=50):
    loo = LeaveOneOut()
    X_arr = np.asarray(X)
    y = np.asarray(y).ravel()
    y_pred = np.zeros(len(y))
    inner_rmse = []
    fold_params = []
    route = 'NestedLOOCV'
    for tr, te in loo.split(X_arr):
        X_tr, X_te = X_arr[tr], X_arr[te]
        y_tr = y[tr]
        nn = len(tr)
        ns = min(5, nn)
        inner = KFold(n_splits=ns, shuffle=True, random_state=42) if ns >= 2 else None
        est = clone(base_estimator)
        if not param_space:
            est.fit(X_tr, y_tr)
            y_pred[te[0]] = float(est.predict(X_te)[0])
            fold_params.append({})
            continue
        if inner is None:
            est.fit(X_tr, y_tr)
            y_pred[te[0]] = float(est.predict(X_te)[0])
            fold_params.append({})
            continue
        n_comb = len(list(ParameterGrid(param_space)))
        if n_comb <= grid_threshold:
            search = GridSearchCV(est, param_space, cv=inner, scoring='neg_mean_squared_error', n_jobs=-1, verbose=0)
            route = f'NestedGrid({n_comb})'
        else:
            ni = min(random_iter, n_comb)
            search = RandomizedSearchCV(
                est, param_distributions=param_space, n_iter=ni, cv=inner,
                scoring='neg_mean_squared_error', random_state=42, n_jobs=-1, verbose=0,
            )
            route = f'NestedRandom({ni}/{n_comb})'
        search.fit(X_tr, y_tr)
        y_pred[te[0]] = float(search.predict(X_te)[0])
        inner_rmse.append(float(np.sqrt(-search.best_score_)))
        fold_params.append(search.best_params_)
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    tune_rmse = float(np.mean(inner_rmse)) if inner_rmse else np.nan
    import json as _j
    from collections import Counter as _C
    keys = [_j.dumps(p, sort_keys=True, default=str) for p in fold_params if p]
    maj = _j.loads(_C(keys).most_common(1)[0][0]) if keys else {}
    final_m = clone(base_estimator)
    if maj:
        final_m.set_params(**maj)
    final_m.fit(X_arr, y)
    return final_m, route, maj, tune_rmse, r2, rmse, y_pred

'''
    if "def nested_loocv_tune_eval" not in s3:
        s3 = s3.replace(
            '# 5. 检查与验证\n',
            NESTED_FN + '\n# 5. 检查与验证\n',
        )
    nb["cells"][3]["source"] = to_nb_source(s3)

    nb["cells"][8]["source"] = to_nb_source(TENSILE_CELL8)
    nb["cells"][9]["source"] = to_nb_source(TENSILE_CELL9)

    s10 = "".join(nb["cells"][10]["source"])
    s10 = s10.replace("an# ============================================", "# ============================================")
    s10 = s10.replace(
        "df_local = pd.read_csv('Final_Features_Ultimate_all_deletion.csv')\n",
        "df_local = df  # 与上文 CSV 读取一致，避免重复路径\n",
    )
    if "from sklearn.model_selection import train_test_split" not in s10:
        s10 = s10.replace(
            "from sklearn.inspection import permutation_importance\n",
            "from sklearn.inspection import permutation_importance\nfrom sklearn.model_selection import train_test_split\n",
        )
    block_fit_perm_old = """    X_df, y = get_clean_local(target_col, feats)
    mname, model = get_best_model_from_upper(target_col)
    model.fit(X_df.values, y)
    print(f"\\n{'='*80}")
    print(f"[{target_name}] 最优算法: {mname} | 样本数: {len(y)} | 特征数: {len(feats)}")
    print(f"{'='*80}")

    # ── 1) Permutation Importance (基于 R2) ──
    print(f"  → 计算 Permutation Importance ...")
    perm_result = permutation_importance(
        model, X_df.values, y, n_repeats=30, random_state=42, scoring='r2'
    )
"""
    block_fit_perm_new = """    X_df, y = get_clean_local(target_col, feats)
    mname, model = get_best_model_from_upper(target_col)
    print(f"\\n{'='*80}")
    print(f"[{target_name}] 最优算法: {mname} | 样本数: {len(y)} | 特征数: {len(feats)}")
    print(f"{'='*80}")

    # ── 1) Permutation Importance（20% 留出集，避免与全训练集重叠）──
    print(f"  → 计算 Permutation Importance ...")
    X_pi_tr, X_pi_val, y_pi_tr, y_pi_val = train_test_split(
        X_df.values, y, test_size=0.2, random_state=42, shuffle=True
    )
    model.fit(X_pi_tr, y_pi_tr)
    perm_result = permutation_importance(
        model, X_pi_val, y_pi_val, n_repeats=30, random_state=42, scoring='r2'
    )
"""
    if block_fit_perm_old in s10:
        s10 = s10.replace(block_fit_perm_old, block_fit_perm_new)

    old_shap_block = """    if mname == 'XGBoost':
        exp = shap.TreeExplainer(model)
        sv = exp.shap_values(X_df)
    elif mname in ('SVR', 'GPR'):
        scaler = model.named_steps['standardscaler']
        inner_key = 'svr' if mname == 'SVR' else 'gaussianprocessregressor'
        inner_model = model.named_steps[inner_key]
        X_sc = scaler.transform(X_df.values)
        bg = shap.kmeans(X_sc, k=30).data
        exp = shap.KernelExplainer(inner_model.predict, bg)
        sv = exp.shap_values(X_sc, nsamples=500)
    else:
        exp = shap.Explainer(model.predict, X_df)
        sv = exp(X_df).values
"""
    new_shap_block = """    if mname == 'XGBoost':
        exp = shap.TreeExplainer(model)
        sv = exp.shap_values(X_pi_val)
    elif mname in ('SVR', 'GPR'):
        scaler = model.named_steps['standardscaler']
        inner_key = 'svr' if mname == 'SVR' else 'gaussianprocessregressor'
        inner_model = model.named_steps[inner_key]
        X_sc = scaler.transform(X_pi_val)
        X_tr_sc = scaler.transform(X_pi_tr)
        bg = shap.kmeans(X_tr_sc, k=min(30, len(X_tr_sc))).data
        exp = shap.KernelExplainer(inner_model.predict, bg)
        sv = exp.shap_values(X_sc, nsamples=min(500, 2000))
    else:
        exp = shap.Explainer(model.predict, pd.DataFrame(X_pi_tr, columns=feats))
        sv = exp(pd.DataFrame(X_pi_val, columns=feats)).values
"""
    if old_shap_block in s10 and "X_pi_val" not in old_shap_block:
        s10 = s10.replace(old_shap_block, new_shap_block)

    nb["cells"][10]["source"] = to_nb_source(s10)

    rep_a12 = (
        "    tuned_model, route_used, best_params, tune_rmse = tune_model_with_dynamic_route(\n"
        "        name, base_models[name], param_spaces[name], X, y\n"
        "    )\n\n"
        "    r2, rmse, _ = evaluate_model(tuned_model, X, y)\n"
    )
    rep_a1314 = (
        "    tuned_model, route_used, best_params, tune_rmse = tune_model_with_dynamic_route(\n"
        "        name, base_models[name], param_spaces[name], X, y\n"
        "    )\n"
        "    r2, rmse, _ = evaluate_model(tuned_model, X, y)\n"
    )
    rep_an = (
        "    tuned_model, route_used, best_params, tune_rmse, r2, rmse, _ = nested_loocv_tune_eval(\n"
        "        name, base_models[name], param_spaces[name], X, y\n"
        "    )\n"
    )
    rep_b = (
        "        tuned_model, route_used, best_params, tune_rmse = tune_model_with_dynamic_route(\n"
        "            name, base_models[name], param_spaces[name], X_cur, y_cur\n"
        "        )\n"
        "        r2, rmse, _ = evaluate_model(tuned_model, X_cur, y_cur)\n"
    )
    rep_bn = (
        "        tuned_model, route_used, best_params, tune_rmse, r2, rmse, _ = nested_loocv_tune_eval(\n"
        "            name, base_models[name], param_spaces[name], X_cur, y_cur\n"
        "        )\n"
    )
    for idx in (12, 13, 14, 20):
        sx = "".join(nb["cells"][idx]["source"])
        sx2 = sx.replace(rep_a12, rep_an) if idx == 12 else sx.replace(rep_a1314, rep_an)
        if sx2 == sx:
            sx2 = sx.replace(rep_b, rep_bn)
        nb["cells"][idx]["source"] = to_nb_source(sx2)

    if "plt.rcParams['font.sans-serif'] = ['SimHei']" in "".join(nb["cells"][1]["source"]):
        s1 = "".join(nb["cells"][1]["source"])
        s1 = s1.replace(
            "plt.rcParams['font.sans-serif'] = ['SimHei'] # 显示中文",
            "plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'Hiragino Sans GB', 'SimHei', 'DejaVu Sans']",
        )
        nb["cells"][1]["source"] = to_nb_source(s1)

    path.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print("Patched", path)


if __name__ == "__main__":
    patch_final(REPO / "forecast_ML_3.8" / "final_3.8.ipynb")
    patch_tensile(REPO / "forecast_ML_1.9" / "alghorithm_multi_1.9_tensile.ipynb")
