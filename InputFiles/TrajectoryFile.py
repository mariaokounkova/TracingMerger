f = open('Trajectories.input', 'w')
f.write('RunBaseDir = Trace;\n')
f.write('Indices = ')
for i in range(100000):
        f.write(str(i)+',')
f.write(';\n')
f.close()

