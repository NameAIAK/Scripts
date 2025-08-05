import json

with open('data.json', 'r') as file:
    data = json.load(file)

successdata=open('successdata.txt', 'a')
successdata.write('dir\tline\tbc\n')

for key, value in data.items():
    dir = key.split('_')[0]
    line=key.split('_')[1]
    bc=key.split('_')[2]

    content = dir+'\t'+line+'\t'+bc+'\n'
    successdata.write(content)