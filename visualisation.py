import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


OBSERV_NUM = 222
COLORS = ('black', 'dimgrey', 'brown', 'red', 'tomato', 'orangered', 'sienna', 'tan', 'orange', 'gold',
          'olive', 'olivedrab', 'yellowgreen', 'forestgreen', 'teal', 'dodgerblue',  'midnightblue',
          'blue', 'slateblue', 'rebeccapurple', 'mediumorchid', 'darkmagenta', 'fuchsia', 'deeppink',
          'palegreen', 'rosybrown', 'goldenrod', 'pink', 'tan')


class Observatory:
    def __init__(self, code, color):
        self.code = code
        self.color = color
        self.times = []
        self.delta_ra = []
        self.delta_dec = []


with open("Data/observatories.txt") as observ_data:
    observ_codes = ['' for _ in range(OBSERV_NUM)]
    observatories = {}
    i = 0
    for row in observ_data.readlines():
        if row[14] != 's':
            code = row[-4:-1]
            if not (code in observ_codes):
                observatories[code] = Observatory(
                    code, COLORS[len(observatories)])
            observ_codes[i] = code
            i += 1

with open("Data/ModelingData.txt") as model_data:
    for i, row in enumerate(model_data.readlines()):
        s_split = row.split()
        code = observ_codes[i]
        t, ra, dec = list(map(float, [s_split[0], s_split[-2], s_split[-1]]))
        observ = observatories[code]
        observ.times.append(t)
        observ.delta_ra.append(ra)
        observ.delta_dec.append(dec)


fig = plt.figure(figsize=(10, 7))

plt.title("Model-Observing delta ra")
plt.xlabel('Time')
plt.ylabel("Delta ra")
plt.grid(linewidth=0.5)
legend_lables = []
for obs in observatories.values():
    plt.scatter(obs.times, obs.delta_ra,
                c=mcolors.CSS4_COLORS[obs.color], label=obs.code)

plt.legend(loc='upper right',  ncol=3)
fig.savefig("graphics/delta_ra.png")

fig.clear()
plt.title("Model-Observing delta dec")
plt.xlabel('Time')
plt.ylabel("Delta dec")
plt.grid(linewidth=0.5)
legend_lables = []
for obs in observatories.values():
    plt.scatter(obs.times, obs.delta_dec,
                c=mcolors.CSS4_COLORS[obs.color], label=obs.code)

plt.legend(loc='upper right', ncol=3)
fig.savefig("graphics/delta_dec.png")
print("Images have been saved in directory './graphics'")