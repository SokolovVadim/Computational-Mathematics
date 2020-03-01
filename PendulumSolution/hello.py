import matplotlib as mpl
import matplotlib.pyplot as plt
import math

dpi = 80
fig = plt.figure(dpi = dpi, figsize = (512 / dpi, 384 / dpi) )
mpl.rcParams.update({'font.size': 10})

plt.axis([0, 300, -1.5, 1.5])

plt.title('Sine & Cosine')
plt.xlabel('x')
plt.ylabel('F(x)')

xs = []
sin_vals = []
cos_vals = []

file = open("output.txt", "r")
data = [map(float, line.split(" ")) for line in file]

x = 0
while x < 300.0:
	print(data.keys[0])
	x += 1
    
    
'''
plt.plot(xs, sin_vals, color = 'blue', linestyle = 'solid',
         label = 'sin(x)')
plt.plot(xs, cos_vals, color = 'red', linestyle = 'dashed',
         label = 'cos(x)')

plt.legend(loc = 'upper right')
fig.savefig('trigan.png')
'''

