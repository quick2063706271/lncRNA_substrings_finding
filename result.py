import pandas as pd

df = pd.read_csv('./result_20210810/results_15.txt', header=None, )

for i in range(16, 26):
    df_n = pd.read_csv('./result_20210810/results_' + str(i) + '.txt', header=None)
    df = df.append(df_n)
df.index = [i for i in range(len(df))]

df.to_csv('./results_ychromo.csv', header=None)
