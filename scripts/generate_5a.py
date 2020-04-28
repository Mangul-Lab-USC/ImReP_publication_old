import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

colors = ['#F7941D', '#00b9f2', '#00a875', '#ecde38', '#0072bc', '#F15a22', '#da6fab']
sns.set_palette(sns.color_palette(colors))
sns.set_context('talk')
sns.set_style('white')

both_df = pd.read_csv('../summary_data/Figure5a_data.csv')
print('df loaded')

both_melt = pd.melt(both_df, id_vars=['sample_pair', 'relationship'], 
                    value_vars=['IGH', 'IGK', 'IGL'], 
                    var_name='chain', value_name='beta_diversity')
print('df melted')

plt.figure(figsize=(10,8))

g=sns.barplot(x="chain", y="beta_diversity", hue='relationship', data = both_melt, edgecolor='black', linewidth=1, facecolor='white', order=['IGH', 'IGK', 'IGL'], hue_order=['same', 'diff'], saturation=1, estimator=np.mean, ci=None, capsize=.1, errcolor='black')
print('outline barplot drawn')

#sns.barplot(x="chain", y="beta_diversity", hue='relationship', data = both_melt, alpha=.25, hue_order=['same', 'diff'], order=['IGH', 'IGK', 'IGL'], saturation=1, estimator=np.mean, ci=95)
#print('color drawn')

sns.stripplot(x="chain", y="beta_diversity", data = both_melt, order=['IGH', 'IGK', 'IGL'], hue_order=['same', 'diff'], jitter=.33, size=4, alpha=.75)
print('points drawn')

plt.ylabel('Beta diversity')

handles, labels = g.get_legend_handles_labels()
plt.legend(handles, ['across tissues within the same individual', 'across tissues from different individuals'], bbox_to_anchor=(0.47, -.4, 0, 0.), frameon=False, ncol=1, handletextpad=0.2, columnspacing=0.8, loc='lower center', prop={'size': 15})
g.xaxis.label.set_size(20)
g.yaxis.label.set_size(20)
g=sns.despine()

plt.savefig("../figures/Figure5_a.pdf", bbox_inches='tight')
plt.savefig("../figures/Figure5_a.png", bbox_inches='tight')
print('plots saved')

