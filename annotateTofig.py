#! usr/bin/env python

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

title = sys.argv[1]

df = pd.read_csv('annotate.log' ,sep='\t')
df = df.drop(12)
df['Log2 Ratio (obs/exp)'] = df['Log2 Ratio (obs/exp)'].astype(float)
print(df)


x = df.loc[:11,['Annotation']]
y = df.loc[:11,['Log2 Ratio (obs/exp)']]

plt.figure(figsize=(10, 7))
plt.bar(x['Annotation'],y['Log2 Ratio (obs/exp)'], width=0.8 ,align='center', alpha=0.5, color='dimgray', capsize=6)
plt.rcParams["font.size"] = 8
plt.title(title)
plt.ylabel('Log2 Ratio (obs/exp)')
plt.xlabel('Annotation')
plt.tight_layout()
plt.savefig(title+'_annotated.png', bbox_inches='tight')