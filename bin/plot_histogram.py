import sys
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource, Legend, Range1d
from bokeh.plotting import figure, save, output_file


def plot_histogram(sim_file, ori_file, outprefix):
    h2_list = []
    output_file(outprefix + ".html")
    with open(sim_file) as source:
        for line in source:
            line = line.rstrip().split()
            h2_list.append(float(line[0]))
    with open(ori_file) as source:
        for line in source:
            line = line.rstrip().split()
            ori_h2 = float(line[0])
    h2_array = np.array(h2_list)
    hist, edges = np.histogram(h2_array, density=False, bins=10)
    p = figure()
    p.output_backend = "svg"
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white")
    p.x_range = Range1d(0, 1)
    p.vspan(x=[ori_h2], line_width=[2], line_color="blue")
    save(p)


if __name__ == "__main__":
    plot_histogram(sys.argv[1], sys.argv[2], sys.argv[3])
