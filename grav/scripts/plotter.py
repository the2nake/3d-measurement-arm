import matplotlib.animation as anim
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

traces = []
with open("out/trace.txt", "r") as file:
    grouping = int(file.readline().strip())
    for _ in range(grouping):
        traces.append([[], []])
    i = 0
    while True:
        line = file.readline()
        if not line:
            break
        split = line.split()
        traces[i % grouping][0].append(float(split[0]))
        traces[i % grouping][1].append(float(split[1]))
        i += 1

lines = []
for i, (x, y) in enumerate(traces):
    color = "red" if i == len(traces) - 1 else "black"
    lines.append(ax.plot(x[0], y[0], color=color, linewidth=0.5)[0])


def update(frame):
    # for each frame, update the data stored on each artist.
    for i, (x, y) in enumerate(traces):
        x_l = x[:frame]
        y_l = y[:frame]
        lines[i].set_xdata(x_l)
        lines[i].set_ydata(y_l)

    return lines


ani = anim.FuncAnimation(fig=fig, func=update, frames=len(traces[0][0]), interval=100)
FFwriter = anim.FFMpegWriter(fps=10)
ani.save("out/animation.mp4", writer=FFwriter)
# plt.show()
