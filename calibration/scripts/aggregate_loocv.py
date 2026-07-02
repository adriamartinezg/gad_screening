import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score, auc, roc_curve
import os
import argparse
import sys
import numpy as np

def apply_extra_gadc_filters(df, tm_min=10, tm_max=12):
    """Structural filters for gadC: Length P10-P90 and TM domains 10-12."""
    if df.empty: return df
    
    if "seq_length" in df.columns:
        positives = df[df["label"] == 1]
        if not positives.empty and len(positives) > 5:
            l_min = positives["seq_length"].quantile(0.10)
            l_max = positives["seq_length"].quantile(0.90)
            df = df[(df["seq_length"] >= l_min) & (df["seq_length"] <= l_max)].copy()

    if "n_tm" in df.columns:
        df = df[(df["n_tm"] >= tm_min) & (df["n_tm"] <= tm_max)].copy()
        
    return df

def load_all_rounds(gene, tsv_list, min_coverage=60.0):
    """Loads and cleans all TSV files from LOOCV rounds."""
    all_dfs = []
    for f in tsv_list:
        if not os.path.exists(f) or os.path.getsize(f) == 0: continue
        try:
            temp_df = pd.read_csv(f, sep="\t")
            if temp_df.empty: continue
            # 'pos' files are 1, 'neg' are 0
            temp_df["label"] = 1 if "_pos" in f else 0
            all_dfs.append(temp_df)
        except Exception as e:
            print(f"Error reading {f}: {e}")

    if not all_dfs: return pd.DataFrame()

    df_combined = pd.concat(all_dfs, ignore_index=True)
    
    score_col = "Score_best_domain" if "Score_best_domain" in df_combined.columns else "score"
    if score_col in df_combined.columns:
        df_combined["bitscore"] = df_combined[score_col]
    else:
        return pd.DataFrame() 
    
    # Coverage filter
    if 'Coverage' in df_combined.columns:
        df_combined = df_combined[df_combined['Coverage'] >= min_coverage].copy()

    # Gene-specific structural filtering
    if gene == "gadC":
        df_combined = apply_extra_gadc_filters(df_combined)
    
    return df_combined

def calc_metrics(df):
    """Calculates TP, FP, Precision, Recall and FPR for all possible thresholds."""
    if df.empty or "label" not in df.columns or len(df["label"].unique()) < 2:
        return pd.DataFrame(columns=["threshold","TP","FP","precision","recall","FPR"])
    
    thresholds = sorted(df["bitscore"].unique(), reverse=True)
    rows = []
    total_p = (df["label"] == 1).sum()
    total_n = (df["label"] == 0).sum()

    for t in thresholds:
        pred = (df["bitscore"] >= t).astype(int)
        tp = ((pred == 1) & (df["label"] == 1)).sum()
        fp = ((pred == 1) & (df["label"] == 0)).sum()
        precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
        recall = tp / total_p if total_p > 0 else 0.0
        fpr = fp / total_n if total_n > 0 else 0.0
        rows.append([t, tp, fp, precision, recall, fpr])

    return pd.DataFrame(rows, columns=["threshold","TP","FP","precision","recall","FPR"])

def select_TC(metrics_df):
    """Heuristic to select the Trusted Cutoff (TC). Priority: Precision >= 0.98."""
    if metrics_df.empty: return 0.0
    mask = (metrics_df.precision >= 0.98) & (metrics_df.recall >= 0.50)
    valid = metrics_df[mask]
    if not valid.empty:
        return valid["threshold"].min()
    return metrics_df["threshold"].max()

def get_metrics_at_threshold(metrics_df, target_t):
    """Retrieves precision/recall for a specific bitscore (closest match)."""
    if metrics_df.empty: return 0.0, 0.0
    idx = (metrics_df['threshold'] - target_t).abs().idxmin()
    row = metrics_df.loc[idx]
    return row['precision'], row['recall']

def plotting(df, gene, output_png):
    """Generates PR and ROC curves."""
    if df.empty or len(df["label"].unique()) < 2:
        plt.figure(); plt.text(0.5, 0.5, "Insufficient data (only one class)", ha='center')
        plt.savefig(output_png); plt.close()
        return 0.0, 0.0

    y_true, y_scores = df["label"], df["bitscore"]
    pr_auc = average_precision_score(y_true, y_scores)
    fpr_vals, tpr_vals, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr_vals, tpr_vals)
    precision_vals, recall_vals, _ = precision_recall_curve(y_true, y_scores)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(recall_vals, precision_vals, color="b", lw=2)
    ax1.set_title(f"{gene} PR Curve (AUC={pr_auc:.3f})")
    ax1.set_xlabel("Recall"); ax1.set_ylabel("Precision"); ax1.grid(True)
    
    ax2.plot(fpr_vals, tpr_vals, color="r", lw=2, label=f"AUC={roc_auc:.3f}")
    ax2.plot([0,1],[0,1],'k--')
    ax2.set_title(f"{gene} ROC Curve")
    ax2.legend(); ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()
    return pr_auc, roc_auc

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--tsvs", nargs="+", required=True)
    parser.add_argument("--output_csv", required=True) 
    parser.add_argument("--output_png", required=True) 
    parser.add_argument("--output_details", required=True) 
    parser.add_argument("--output_threshold_plot", required=True) 
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    # 1. Load Data
    df = load_all_rounds(args.gene, args.tsvs)

    if df.empty:
        print(f"WARNING: No data for {args.gene}. Generating empty files.")
        pd.DataFrame().to_csv(args.output_csv, index=False)
        sys.exit(0)

    # 2. Compute Threshold Metrics
    metrics_df = calc_metrics(df)
    if args.output_details:
        metrics_df.to_csv(args.output_details, sep="\t", index=False)

    # 3. Select TC and robustness (TC +/- 2)
    tc = select_TC(metrics_df)
    prec_tc, rec_tc = get_metrics_at_threshold(metrics_df, tc)
    prec_m2, rec_m2 = get_metrics_at_threshold(metrics_df, tc - 2.0)
    prec_p2, rec_p2 = get_metrics_at_threshold(metrics_df, tc + 2.0)

    # 4. Plotting
    if args.output_threshold_plot:
        plt.figure(figsize=(8, 6))
        plt.plot(metrics_df["threshold"], metrics_df["precision"], label="Precision", color="green", lw=2)
        plt.plot(metrics_df["threshold"], metrics_df["recall"], label="Recall", color="orange", lw=2)
        plt.axvline(x=tc, color="red", linestyle="--", label=f"TC: {tc:.1f}")
        plt.axvline(x=tc-2, color="gray", linestyle=":", alpha=0.5, label="TC-2")
        plt.axvline(x=tc+2, color="gray", linestyle=":", alpha=0.5, label="TC+2")
        plt.title(f"Threshold Stability: {args.gene}")
        plt.xlabel("Bitscore"); plt.ylabel("Score"); plt.legend(); plt.grid(alpha=0.3)
        plt.savefig(args.output_threshold_plot, dpi=300); plt.close()

    pr_auc, roc_auc = plotting(df, args.gene, args.output_png)

    # 5. Final Summary CSV
    res_df = pd.DataFrame([{
        "gene": args.gene, 
        "TC": tc, 
        "PR_AUC": pr_auc, 
        "ROC_AUC": roc_auc, 
        "n_pos": (df["label"] == 1).sum(), 
        "n_neg": (df["label"] == 0).sum(),
        "Prec@TC": round(prec_tc, 4),
        "Rec@TC": round(rec_tc, 4),
        "Prec@TC-2": round(prec_m2, 4),
        "Rec@TC-2": round(rec_m2, 4),
        "Prec@TC+2": round(prec_p2, 4),
        "Rec@TC+2": round(rec_p2, 4)
    }])
    
    res_df.to_csv(args.output_csv, index=False)
    print(f"Summary for {args.gene} saved.")
