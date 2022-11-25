#Define results directories and style sheets
#Define colors and plot helper functions
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import os
#plt.style.use('../scripts/nmjact.mplstyle')
#plt.style.use('scripts/nmjact.mplstyle')
indir = '../'
plt.style.use(os.path.join(indir, 'scripts/nmjact.mplstyle'))

#config file copied from snakemake will be here:
config_file = os.path.join(indir, 'scripts/config.yml')
with open(config_file, 'r') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

results_dir = config['results_dir']
results_dir_inctest = config['results_dir_inctest']
results_dir_pherson = config['results_dir_pherson']
geo_outdir = config['geo_outdir']
pherson_dir = config['pherson_dir']
inspect_dir = config['inspect_dir']
gffutils_db = config['gffutils_db']
out_fmt = config['out_fmt']
out_dpi = config['out_dpi']

if not os.path.exists(geo_outdir):
    os.makedirs(geo_outdir)

carto12_hex = ['#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310',
'#008695', '#CF1C90', '#f97b72', '#4b4b8f', '#A5AA99']

#https://gist.github.com/JoachimGoedhart/b2f91652de2b9e3b393c6c28be843e00
#these are the tol_muted colors
tol_muted = {'indigo': '#332288', 'cyan':'#88CCEE', 'teal':'#44AA99', 'green':'#117733',
'olive':'#999933', 'sand':'#DDCC77', 'rose':'#CC6677', 'wine':'#882255', 'purple':'#AA4499',
 'pale_grey':'#DDDDDD', 'grey':'#BBBBBB'}

tol_bright = {'red':'#EE6677', 'green':'#228833', 'blue':'#4477AA', 'yellow':'#CCBB44', 'cyan':'#66CCEE', 'purple':'#AA3377', 'grey':'#BBBBBB'}
color_dict = tol_bright.copy()
#could maybe swap red and cyan
ordered_colors = [tol_bright[i] for i in ['grey', 'purple', 'blue', 'red', 'green', 'cyan', 'yellow']]
#tol_muted = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499', '#DDDDDD']

##selected_colors = carto12_hex
selected_colors = ordered_colors

# alpha = 0.3
# light_colors = []
# for c in ordered_colors:
#     rgba = [int(c.lstrip('#')[i:i+2], 16)/255 for i in (0,2,4)]
#     rgba.append(alpha)
#     light_colors.append(rgba)

#put the grey color first
#selected_colors.reverse()
#What happens if we do this?
sns.set_palette(sns.color_palette(selected_colors))
#Does this set prop cycle?
#mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=ordered_colors)

#column widths, converting to inches for matplotlib
#these are from the Cell family journals for now:
# Gutter, to add in mm between the subplots
gutter = 1/25.4

scol = 85/25.4
mcol = 114/25.4
dcol = 174/25.4
# To account for the gutter space, multiply gutter*(n-1)
sfig = (scol - gutter*1)/2
dfig = (dcol - gutter*3)/4
#I think generally the subpanel will be a quarter of page width, so
#single_col/2 or double_col/4 or mcol/3 (1.5)

figsizes = {'single_col': (sfig, sfig), 'double_col': (dfig, dfig)}
#mpl.rcParams["legend.edgecolor"] = 'k'
# https://stackoverflow.com/questions/38359008/matplotlib-is-there-an-rcparams-way-to-adjust-legend-border-width
# https://stackoverflow.com/questions/27698377/how-do-i-make-sans-serif-superscript-or-subscript-text-in-matplotlib
mpl.rcParams['patch.linewidth'] = 0.75
mpl.rcParams['mathtext.default'] = 'regular'
