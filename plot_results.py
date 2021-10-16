import os
import sys
from pandas import read_table
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    print("Usage: %s <input file>" % sys.argv[0])

table = read_table(input_file)
print(table)
j = table["j"]
var = table["VAR"]
mag = table["Mw"]
qty = table["QUALITY"]
pdc = table["DC"]
pclvd = table["CLVD"]
piso = table["ISO"]

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(j.values, var.values, color='red', label='Variance')
ax2.plot(j.values, qty.values, color='blue', label='Quality')

ax1.yaxis.set_major_locator(MultipleLocator(20))
ax1.yaxis.set_minor_locator(MultipleLocator(10))
ax1.xaxis.set_major_locator(MultipleLocator(2))
ax1.xaxis.set_minor_locator(MultipleLocator(1))
ax2.yaxis.set_major_locator(MultipleLocator(2))
ax2.yaxis.set_minor_locator(MultipleLocator(1))
ax1.set_xlabel("Scale")
ax2.set_ylabel("Quality")
ax1.set_ylabel("Varianc %")
ax1.set_ylim([-10, 110])
ax2.set_ylim([-1, 5])


ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

left, bottom, width, height = [0.17, 0.65, 0.25, 0.25]
ax3 = fig.add_axes([left, bottom, width, height])
ax3.plot(j.values, 1-pclvd.values/pdc.values, color='green', label='DC ratio')
# ax3.set_ylim([-.2, 1.2])
ax3.legend(loc="lower right", fontsize=6)
ax3.xaxis.set_major_locator(MultipleLocator(2))
ax3.xaxis.set_minor_locator(MultipleLocator(1))
ax3.tick_params(axis='both', which='both', labelsize=6)
ax3.set_ylabel("DC ratio", fontsize=6)
ax3.set_xlabel("Scales", fontsize=6)
#ax3.set_facecolor("none")

left, bottom, width, height = [0.4, 0.2, 0.25, 0.25]
ax4 = fig.add_axes([left, bottom, width, height])
ax4.plot(j.values, mag.values, color='purple', label='Mw')
# ax3.set_ylim([-.2, 1.2])
ax4.legend(loc="upper right", fontsize=6)
ax4.xaxis.set_major_locator(MultipleLocator(2))
ax4.xaxis.set_minor_locator(MultipleLocator(1))
ax4.yaxis.set_major_locator(MultipleLocator(1))
ax4.yaxis.set_minor_locator(MultipleLocator(.5))
ax4.set_ylim([0, 10])
ax4.set_ylabel("Magnitude Mw", fontsize=6)
ax4.set_xlabel("Scales", fontsize=6)
ax4.tick_params(axis='both', which='both', labelsize=6)

plt.tight_layout()
plt.savefig(os.path.join(".\output\synthetic_test.eps"))
plt.show()
