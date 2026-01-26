import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_path = './temp_2_processed.xlsx'
sheet_name = 'Sheet2'
df = pd.read_excel(file_path, sheet_name=sheet_name)

print(df.head())

df = df[['organism', 'count']]

df.dropna(subset=['organism', 'count'], inplace=True)

organisms = df['organism']
counts = df['count']

x_pos = np.arange(len(organisms)) * 3  

plt.figure(figsize=(12, 6))  
plt.bar(x_pos, counts, color='#cc0000')  

plt.title('Number of Selenoproteins', fontsize=16)
plt.xlabel('Organism', fontsize=12)
plt.ylabel('Count', fontsize=12)

plt.tick_params(axis='y', labelsize=14) 

plt.xticks([]) 

plt.ylim(0, max(counts) + 10)

y_ticks = np.arange(0, max(counts) + 10, 5) 
plt.yticks(y_ticks)

plt.grid(axis='y', linestyle='--', alpha=0.7)  

plt.tight_layout()

plt.savefig('./bar_chart_sort_no_clade_new.pdf', format='pdf', dpi=300)
plt.show()