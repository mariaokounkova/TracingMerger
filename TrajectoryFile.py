f = open('Trajectories.input', 'w')
f.write('RunBaseDir = Kerr;\n')
f.write('Indices = ')
for i in range(10000):
        f.write(str(i)+',')
f.write(';\n')
f.close()

