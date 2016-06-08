import pandas as pd

sec_file = 'uniprot_sec_ac.txt'
lines = open(sec_file, 'rt').readlines()

for i, l in enumerate(lines):
   if l.startswith('Secondary AC'):
       entry_lines = lines[i+2:]

sec_id = []
prim_id = []

for l in entry_lines:
   s, p = l.split()
   sec_id.append(s)
   prim_id.append(p)

d = {'Secondary_ID': sec_id, 'Primary_ID': prim_id}

df = pd.DataFrame(data=d)

df.to_csv('Uniprot_sec_to_prim.csv', sep='\t', index=False) 
