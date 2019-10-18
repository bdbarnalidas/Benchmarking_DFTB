import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gc

x_axis = [11, 20, 41, 80, 167, 314]
y_axis1 = [1.34, 5.43, 6.21, 12.28, 30.90, 92.50]  # DFTB+
y_axis2 = [46.90, 91.92, 124.72, 160.50, 215.48, 302.93]  # Python interface

plt.plot(x_axis, y_axis1, color='blue')
plt.plot(x_axis, y_axis2, color='red')

plt.grid()
blue_patch = mpatches.Patch(color='blue', label='DFTB+')
red_patch = mpatches.Patch(color='red', label='Python interface')

plt.figlegend(handles=[blue_patch, red_patch], loc='upper right', fontsize=9)
# plt.title('Benchmark_graph', fontsize=20)
plt.xlabel('#Atoms', fontsize=16)
plt.ylabel('Time (in seconds)', fontsize=16)
# plt.show()
plt.savefig('Benchmark_graph', bbox_inches='tight', dpi=100)
plt.close()
gc.collect()



