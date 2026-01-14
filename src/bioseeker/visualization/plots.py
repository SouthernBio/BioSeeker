import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_conservation_regression(bicodon_df: pd.DataFrame, slope: float, r_squared: float, output_path: str):
    """
    Plots Observed vs Expected conservation rates and the regression line.
    """
    plt.figure(figsize=(10, 8))
    
    x = bicodon_df['ExpectedRate']
    y = bicodon_df['ConservationRate']
    
    plt.scatter(x, y, alpha=0.5, label='Bicodons')
    
    # Plot regression line
    max_val = max(x.max(), y.max())
    line_x = np.linspace(0, max_val, 100)
    line_y = slope * line_x
    
    plt.plot(line_x, line_y, color='red', label=f'Fit (m={slope:.2f}, $R^2$={r_squared:.2f})')
    
    # Plot identity line (y=x) for reference
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Expected=Observed')
    
    plt.xlabel('Expected Conservation Rate')
    plt.ylabel('Observed Conservation Rate')
    plt.title('Bicodon Conservation: Observed vs Expected')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(output_path)
    plt.close()
