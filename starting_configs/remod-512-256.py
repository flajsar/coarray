i=1
f = open("256.txt", "w")
array=[]
warray=[]
with open('starting_config_TIP4PEW-512.txt') as file:
    for line in file.readlines():
      array.append(line)

for mol in array:
    if i%8<=4 and i%8>0:
        warray.append(mol)
    i+=1
for i in warray:
    f.write(i)
                
        
