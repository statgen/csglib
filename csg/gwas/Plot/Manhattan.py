import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd


def compress_axis(points, diff):
   points_it = iter(sorted(points))
   prev_point = next(points_it)
   yield prev_point
   for point in points_it:
      if point - prev_point > diff:
         prev_point = point
         yield prev_point


def compress_manhattan(data, chromosome_bounds, dpi, width_inch, height_inch, y_break):
   width_px = width_inch * dpi
   height_px = height_inch * dpi
   min_max = data['pvalue'].aggregate([np.min, np.max])
   x_per_px = sum(chromosome_bounds['amax'] - chromosome_bounds['amin']) / width_px
   if y_break:
      y_per_px = min((y_break[0] - min_max['amin']) / (height_px * 0.75), (min_max['amax'] - y_break[1]) / (height_px * 0.25))
   else:
      y_per_px = (min_max['amax'] - min_max['amin']) / height_px

   chrom = None
   center_bp = None
   end_bp = None
   pvalues = []
   compressed = []
   for row in data.itertuples():
      if row.chrom != chrom or row.position > end_bp:
         if pvalues:
            for pvalue in compress_axis(pvalues, y_per_px):
               compressed.append([chrom, center_bp, pvalue])
            pvalues = []
         chrom = row.chrom
         end_bp = row.position + x_per_px
         center_bp = row.position + x_per_px / 2
      pvalues.append(row.pvalue)
   if pvalues:
      for pvalue in compress_axis(pvalues, y_per_px):
         compressed.append([chrom, center_bp, pvalue])
      pvalues = []
   return pd.DataFrame(compressed, columns = ('chrom', 'position', 'pvalue'))


def manhattan_plotter_ybreak(figure, title, data, chromosome_bounds, y_break):
   gs = gridspec.GridSpec(2, 1, height_ratios = [1, 3])
   ax1 = figure.add_subplot(gs[0, 0])
   ax2 = figure.add_subplot(gs[1, 0])

   figure.text(0.09, 0.5, r'$-\log_{10}$(p-value)', va='center', rotation='vertical')
   ax2.set_xlabel('Chromosomes')

   if title:
      ax1.set_title(title)

   chromosomes = sorted(data['chrom'].unique(), key = lambda x: float(x) if x.isdigit() else float('inf'))
   x_per_px = sum(chromosome_bounds['amax'] - chromosome_bounds['amin']) / (figure.get_size_inches()[0] * figure.get_dpi())

   start_y = 0
   ticks = []
   for name in chromosomes:
      bounds = chromosome_bounds[chromosome_bounds['chrom'] == name]
      if bounds.empty:
         continue
      start = bounds['amin'].values[0]
      end = bounds['amax'].values[0]
      ax1.scatter(start_y + data[data['chrom'] == name]['position'], data[data['chrom'] == name]['pvalue'], s = 3)
      ax2.scatter(start_y + data[data['chrom'] == name]['position'], data[data['chrom'] == name]['pvalue'], s = 3)
      ticks.append(start_y + (end - start) / 2)
      start_y += end + 3 * x_per_px

   ax1.set_ylim(y_break[1], ax1.get_ylim()[1])
   ax1.set_xlim(0, start_y)

   ax2.set_ylim(0, y_break[0])
   ax2.set_xlim(0, start_y)

   if ax1.get_yticks()[0] == y_break[1]:
      ax1.set_yticks(ax1.get_yticks()[1:])

   if ax2.get_yticks()[-1] == y_break[0]:
      ax2.set_yticks(ax2.get_yticks()[:-1])

   ax1.spines['bottom'].set_visible(False)
   ax1.set_xticks([])
   ax2.spines['top'].set_visible(False)
   ax2.set_xticks(ticks)
   ax2.xaxis.set_ticks_position('bottom')
   ax2.set_xticklabels([str(i) for i in xrange(1, 23)])

   d = .005  # how big to make the diagonal lines in axes coordinates
   # arguments to pass to plot, just so we don't keep repeating them
   kwargs = dict(transform = ax1.transAxes, color = 'k', clip_on = False)
   ax1.plot((-d, +d), (-d - 0.01, +d + 0.01), **kwargs)        # top-left diagonal
   ax1.plot((1 - d, 1 + d), (-d - 0.01, +d + 0.01), **kwargs)  # top-right diagonal
   kwargs.update(transform = ax2.transAxes)  # switch to the bottom axes
   ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
   ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
   ax2.axhline(y = -np.log10(5 * 10**-8), color = 'red', linestyle = '--', linewidth = 1)
   plt.subplots_adjust(hspace = 0.05)


def manhattan_plotter(figure, title, data, chromosome_bounds):
   ax = figure.add_subplot(1, 1, 1)

   ax.set_ylabel(r'$-\log_{10}$(p-value)')
   ax.set_xlabel('Chromosomes')

   if title:
      ax.set_title(title)

   chromosomes = sorted(data['chrom'].unique(), key = lambda x: float(x) if x.isdigit() else float('inf'))
   x_per_px = sum(chromosome_bounds['amax'] - chromosome_bounds['amin']) / (figure.get_size_inches()[0] * figure.get_dpi())

   start_y = 0
   ticks = []
   for name in chromosomes:
      bounds = chromosome_bounds[chromosome_bounds['chrom'] == name]
      if bounds.empty:
         continue
      start = bounds['amin'].values[0]
      end = bounds['amax'].values[0]
      ax.scatter(start_y + data[data['chrom'] == name]['position'], data[data['chrom'] == name]['pvalue'], s = 3)
      ticks.append(start_y + (end - start) / 2)
      start_y += end + 3 * x_per_px

   ax.set_xlim(0, start_y)
   ax.set_ylim(0, ax.get_ylim()[1])
   ax.set_xticks(ticks)
   ax.set_xticklabels(chromosomes)

   ax.axhline(y = -np.log10(5 * 10**-8), color = 'red', linestyle = '--', linewidth = 1)


def create_manhattan(data, file_name, y_break, title, dpi, width_inches, height_inches):
   data.sort_values(['chrom', 'position'], axis = 0, ascending = [True, True], inplace = True)
   data.reset_index(drop = True, inplace = True)
   data['pvalue'] = -np.log10(data['pvalue'])
   chromosome_bounds = data.groupby(['chrom'], sort = False)['position'].aggregate([np.min, np.max]).reset_index()

   matplotlib.style.use('seaborn-paper')

   fig = plt.figure(figsize = (width_inches, height_inches), dpi = dpi)


   compressed = compress_manhattan(data, chromosome_bounds, dpi, width_inches, height_inches, y_break)
   #print len(compressed.index), len(data.index)
   if y_break:
      manhattan_plotter_ybreak(fig, title, compressed, chromosome_bounds, y_break)
   else:
      manhattan_plotter(fig, title, compressed, chromosome_bounds)

   plt.savefig(file_name, dpi = dpi, format = "png", bbox_inches = 'tight')
   plt.show()
