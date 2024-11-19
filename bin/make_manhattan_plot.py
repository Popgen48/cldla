import sys
import os.path
import pandas as pd
import yaml
from yaml.loader import SafeLoader
import math
from bokeh.models import (
    HoverTool,
    ColumnDataSource,
    Span,
    CustomJS,
    TapTool,
)
from bokeh.plotting import figure, save, output_file
from bokeh.models.tickers import FixedTicker


class InteractiveManhattanplot:
    def __init__(self, chrom_cord_val, yml_file, chrom_cutoff_f, window_size, y_label, outprefix):
        self.outprefix = outprefix  # prefix of the output file
        self.chrom_cord_val = chrom_cord_val  # input file should be sorted by chromosome id i.e. 1,2,3,4,...,n
        self.yml_file = yml_file
        self.window_size = int(window_size)  # use for ensembl link
        # self.min_score = 0  # important to set the min cordi of y axis
        # self.max_score = 0  # important to set the maxi cordi of y axis
        self.y_label = y_label  # label will determine the hovering points, has to be one of the string present in self.upper_hover or tajimas_d
        self.upper_hover = ["fst_values", "LRT_values", "-log10_P"]
        self.axis_chrom_dict = (
            {}
        )  # dictionary with greatest coordinates of a respective chromosome as key and chromosome its value --> use to replace major X-axis with chromosome name
        self.chrom_cutoff_f = chrom_cutoff_f
        self.chrom_cutoff_d = {}
        self.line_plot_d = {}

    def create_chrom_cutoff_dict(self):
        lrt_l = []
        if os.path.isfile(self.chrom_cutoff_f):
            with open(self.chrom_cutoff_f) as source:
                for line in source:
                    line = line.rstrip().split()
                    self.chrom_cutoff_d[line[0]] = line[1]
                    lrt_l.append(float(line[1]))
        self.cutoff = max(lrt_l)

    def read_yml_file(self):
        """
        read yml file and set the parameters corresponding to the interactive plot
        """
        with open(self.yml_file, "r") as p:
            params = yaml.load(p, Loader=SafeLoader)
        self.figure_width = params["width"]  # determine the width of the plot
        self.figure_height = params["height"]  # determine the height of the plot
        self.chrom_label_orientation = params[
            "chrom_label_orientation"
        ]  # determine the angle at which the label will be placed below the X-axis
        self.legend_font_size = params["legend_font_size"]  # determine the font size of the legend
        self.label_font_size = params["label_font_size"]  # determine the font size of the label
        self.fil_alpha = params["fil_alpha"]  # determine the intensity of the color filled in the circle
        self.ensembl_link = params["ensembl_link"]  # source link of the ensembl reference
        self.color_list = params["color_list"]  # color list

    def add_nonsigni_df(self, pd1):
        """
        plot the data of non-significant windows
        """
        source = ColumnDataSource(pd1)
        s = self.p.circle(
            "cum_cord",
            "p_val",
            color="col",
            fill_alpha=self.fil_alpha,
            source=source,
        )

    def add_signi_df(self, pd2):
        """
        plot the data of significant windows
        """
        source_t = ColumnDataSource(pd2)
        t = self.p.circle(
            "cum_cord",
            "p_val",
            color="col",
            fill_alpha=self.fil_alpha,
            source=source_t,
        )
        if self.ensembl_link != "none":
            toolt = """
                <div>
                    <a href=@link{safe}}>@link{safe}</a>
                </div>
                """
            hvr = HoverTool(renderers=[t], tooltips=toolt)
            self.p.add_tools(hvr)

            tap_cb = CustomJS(
                code="""
                              var l = cb_data.source.data['link'][cb_data.source.inspected.indices[0]]
                              window.open(l)
                              """
            )
            tapt = TapTool(renderers=[t], callback=tap_cb, behavior="inspect")
            self.p.add_tools(tapt)
        else:
            hover = HoverTool(tooltips=[("chrom", "@chrom"), ("cord", "@cord")], renderers=[t])
            self.p.add_tools(hover)
            self.p.add_tools(HoverTool(tooltips=[("chrom", "@chrom"), ("cord", "@cord")]))

    def format_plot(self):
        """
        format the plot using the parameters set in yml file
        """
        self.p.xaxis.major_label_text_font_size = self.label_font_size
        self.p.xaxis.ticker = FixedTicker(ticks=list(self.axis_chrom_dict.keys()))
        self.p.xaxis.major_label_overrides = {k: v for k, v in self.axis_chrom_dict.items()}
        self.p.xaxis.major_tick_out = 15
        self.p.xaxis.major_label_orientation = math.radians(self.chrom_label_orientation)
        hline = Span(
            location=self.cutoff,
            dimension="width",
            line_color="red",
            line_width=3,
        )
        self.p.renderers.extend([hline])
        self.p.vspan(x=list(self.axis_chrom_dict.keys()), line_color="black")
        self.p.grid.visible = False
        self.p.xaxis.axis_label = "chromosome_id"
        self.p.yaxis.axis_label = self.y_label

    def sel_out_to_list(self):
        df_list = []  # store the record for which the ensembl database will not be linked
        df_list_t = []  # store the record for which the ensembl database will be linked
        chrom_list = []  # store the chromosome id
        with open(self.chrom_cord_val) as source:
            for line in source:
                line = line.rstrip().split()
                if (
                    line[-1] != "-nan"
                ):  # this takes into account "nan" values generated by the vcftools output for tajimas_d, fst and pi
                    tmp_list = []
                    if line[0] not in chrom_list:
                        if len(chrom_list) != 0:
                            self.axis_chrom_dict[cum_cord_1] = chrom_list[
                                -1
                            ]  # set starting cord of new chromosome as last value of the previous chrom, for first chrom?
                        chrom_list.append(line[0])
                        color = self.color_list[0] if len(chrom_list) % 2 == 0 else self.color_list[1]
                    cum_cord_1 = (
                        int(list(self.axis_chrom_dict.keys())[-1])
                        + int(
                            line[1]
                        )  # to maintain continuity, add the new chrom's cordinates to the last cord of previous chrom
                        if len(chrom_list) > 1
                        else int(line[1])  # if its first chromosome, keep the cord as is
                    )
                    chrom = line[0]
                    cord = line[1]
                    p_val = float(line[-1])
                    tmp_list = [
                        chrom,
                        cord,
                        cum_cord_1,
                        p_val,
                        color,
                    ]  # original cord ('cord') -> link record to ensembl, cum_cord_1 is updated cord, line 143
                    if self.y_label not in self.upper_hover:
                        if p_val < float(self.cutoff) and self.ensembl_link != "none":
                            tmp_list.append(
                                self.ensembl_link + chrom + ":" + cord + "-" + str(int(cord) + self.window_size)
                            )
                        df_list.append(tmp_list[:]) if p_val > float(self.cutoff) else df_list_t.append(tmp_list[:])
                    else:
                        if p_val > float(self.cutoff) and self.ensembl_link != "none":
                            tmp_list.append(
                                self.ensembl_link + chrom + ":" + cord + "-" + str(int(cord) + self.window_size)
                            )
                        df_list_t.append(tmp_list[:]) if p_val > float(self.cutoff) else df_list.append(tmp_list[:])
            self.axis_chrom_dict[cum_cord_1] = chrom_list[-1]
        pd1 = pd.DataFrame(df_list, columns=["chrom", "cord", "cum_cord", "p_val", "col"])
        pd2 = pd.DataFrame(
            df_list_t,
            columns=(
                ["chrom", "cord", "cum_cord", "p_val", "col", "link"]
                if self.ensembl_link != "none"
                else ["chrom", "cord", "cum_cord", "p_val", "col"]
            ),
        )
        return pd1, pd2

    def prepare_lineplot_dict(self):
        if len(self.chrom_cutoff_d) > 0:
            cnt = 0
            for k, v in self.axis_chrom_dict.items():
                if cnt == 0:
                    self.line_plot_d[1] = self.chrom_cutoff_d[v]
                    cnt += 1
                else:
                    self.line_plot_d[previous_pos + 1] = self.chrom_cutoff_d[v]
                self.line_plot_d[k] = self.chrom_cutoff_d[v]
                previous_pos = k

    def draw_line_plot(self):
        if len(self.chrom_cutoff_d) > 0:
            self.p.line(list(self.line_plot_d.keys()), list(self.line_plot_d.values()), line_width=4)

    def main_func(self):
        self.read_yml_file()
        self.p = figure(width=self.figure_width, height=self.figure_height)
        self.p.output_backend = "svg"
        output_file(self.outprefix + ".html")
        self.create_chrom_cutoff_dict()
        pd1, pd2 = self.sel_out_to_list()
        self.add_nonsigni_df(pd1)
        self.add_signi_df(pd2)
        self.prepare_lineplot_dict()
        self.draw_line_plot()
        self.format_plot()
        save(self.p)


if __name__ == "__main__":
    obk = InteractiveManhattanplot(
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
    )
    obk.main_func()
