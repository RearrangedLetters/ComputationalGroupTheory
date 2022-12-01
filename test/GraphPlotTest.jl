using GraphPlot
using Compose

using Graphs
h = watts_strogatz(50, 6, 0.3)
gplot(h)