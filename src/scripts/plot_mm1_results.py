import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Configure plot style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['legend.fontsize'] = 12

def load_data(data_dir: Path) -> pd.DataFrame:
    csv_path = data_dir / 'mm1_varying_lambda.csv'
    df = pd.read_csv(csv_path)
    return df

def plot_theoretical_curve():
    _, ax = plt.subplots(figsize=(10, 6))

    rho = np.linspace(0.01, 0.99, 500)
    L = rho / (1 - rho)

    ax.plot(rho, L, 'b-', linewidth=2.5, label='L = ρ/(1-ρ)')

    ax.axvline(x=1.0, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
               label='Stability boundary (ρ = 1)')

    key_rhos = [0.5, 0.7, 0.9, 0.95]
    for r in key_rhos:
        L_val = r / (1 - r)
        ax.plot(r, L_val, 'ro', markersize=8)
        ax.annotate(f'ρ={r}, L={L_val:.1f}', xy=(r, L_val),
                   xytext=(r-0.15, L_val+2), fontsize=10,
                   arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5))

    ax.set_xlabel('Traffic Intensity (ρ = λ/μ)', fontsize=14)
    ax.set_ylabel('Average Queue Length (L)', fontsize=14)
    ax.set_title('M/M/1 Queue: Theoretical Average Queue Length', fontsize=16)
    ax.legend(loc='upper left', fontsize=12)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 25)
    ax.grid(True, alpha=0.3)

    ax.text(0.05, 20, r'$L = \frac{\rho}{1-\rho} = \frac{\lambda}{\mu - \lambda}$',
            fontsize=16, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.show()

def print_summary_statistics(df: pd.DataFrame):
    print("\n" + "="*60)
    print("M/M/1 Queue Simulation Results Summary")
    print("="*60)

    print(f"\nNumber of experiments: {len(df)}")
    print(f"Rho range: {df['Rho'].min():.2f} to {df['Rho'].max():.2f}")

    stable_df = df[df['Rho'] <= 0.9]
    print(f"\n--- Stable Region (ρ ≤ 0.9) ---")
    print(f"Mean relative error: {stable_df['Relative_Error'].mean():.2f}%")
    print(f"Max relative error: {stable_df['Relative_Error'].max():.2f}%")
    print(f"Experiments within 5% error: {(stable_df['Relative_Error'] <= 5).sum()}/{len(stable_df)}")

    high_df = df[df['Rho'] > 0.9]
    print(f"\n--- High Utilization Region (ρ > 0.9) ---")
    print(f"Mean relative error: {high_df['Relative_Error'].mean():.2f}%")
    print(f"Max relative error: {high_df['Relative_Error'].max():.2f}%")
    print("Note: Higher errors expected due to finite simulation time")

    print("\n" + "="*60)

def main():
    # Paths
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / 'data'
    output_dir = script_dir.parent / 'plots'

    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)

    # Load data
    csv_path = data_dir
    print(f"Loading data from {csv_path}")
    df = load_data(csv_path)

    # Print summary
    print_summary_statistics(df)

    # Generate plots
    plot_theoretical_curve()

    print(f"\nAll plots saved to {output_dir}")

if __name__ == '__main__':
    main()
