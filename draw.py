import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
import glob
import os
import pandas as pd
from energydiagram import ED

diagram = ED("")
diagram.dimension = 12 	# length of horizontal lines
diagram.space = 10		# distance between levels
diagram.aspect = 0.7

level_width = 3.5
link_width = 2.5

dft_alpha = 1

title = 'Gibbs Free Energy [kcal/mol]'

database_name = 'energies.xlsx'
column_to_plot = 'ΔG, kcal/mol'
output_name = f'out/dlpmo-b_s2o'

name_labels = ['Reagents', 'A', 'TS1', 'B', 'TS2', 'C', 'Products']
colors = ['#49B6FF', '#FF6563', '#FCBA04', '#38E3E8', '#77F21D', '#6530F2', 'maroon', 'gold', 'navy']
line_styles = ['solid', 'dashed', 'dashdot', 'dotted']

show_plot = True
export_svg = True
export_png = False
export_pdf = True

draw_images = True
image_dir = f'mols/bs2o/*.png'
image_pos = [15, -30, 22, -25, 18, -25, 20]
image_size = [0.6, 1, 1, 1, 1, 1, 0.6]

draw_e_box = True
e_box_background = 'snow'
e_box_width = 6
e_box_textsize = 16
e_box_offset = 3.5

legend_pos = 0

draw_labels = True
label_textsize = 21
label_offset = 8

squeeze_offset = 25

plt.rcParams["figure.figsize"] = [19.2,10.8]
plt.rcParams['font.family'] = 'CMU Serif'
plt.rcParams['font.serif'] = 'CMU Serif'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'CMU Serif:roman'
plt.rcParams['mathtext.it'] = 'CMU Serif:italic'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24

if draw_images:
	images = []
	for file in glob.glob(image_dir):
		images.append(mpimg.imread(file))

directory_path = os.path.dirname(output_name)
if not os.path.exists(directory_path):
	try:
		os.makedirs(directory_path)
	except OSError as e:
		pass

df = pd.read_excel(database_name)
df = df.drop(df.index[11::12], 0).drop(['Datasource', 'Filenames'], axis=1)

for i, c in enumerate(df.columns):
	if i==0 or i>=7: df[c].fillna(method='ffill', inplace=True)

df = df[(df['Plot']=='Yes') & (df['Method']=='DLPNO//r^2SCAN-3c')].iloc[:,:-1]
df = df.set_index(['Index']+df.columns[8:].tolist()+['Molecule'])

labels = []
c_i = -1
prev_ind = -1
for i, (k, v) in enumerate(df.groupby(level=[0,1,2,3,4])):
	ind, cat, met, sol, add = k

	energies = v[column_to_plot].iloc[3:10]
	labels.append(r'$\rm{%s %s}$' % (met.replace("-", u"\u2212").replace("w", r"\it{\omega}\rm").replace(" ", r"\ "), r'\ (%s)' % sol if sol != 'No' else ''))

	if ind != prev_ind:
		c_i += 1
		prev_ind = ind

	for j, e in enumerate(energies):
		if pd.isna(e):
			continue

		diagram.add_level(e, "", position=j, top_text="",
						  color=colors[i%len(colors)], alpha=1 if 'DLPNO' in met else dft_alpha)
		if j > 0:
			c_id = len(diagram.energies) - 1
			diagram.add_link(c_id - 1, c_id, colors[i%len(colors)], line_styles[i%len(line_styles)], linewidth=link_width)

diagram.add_labels(labels)
diagram.plot(level_width=level_width)

diagram.ax.set_ylabel(title)
if legend_pos == 0:
	diagram.ax.legend(loc='upper right', bbox_to_anchor=(1.05, 1.05))
elif legend_pos == 1:
	diagram.ax.legend(loc='lower left')

l_x, u_x = diagram.ax.get_xlim()
l_y, u_y = diagram.ax.get_ylim()
diagram.ax.set_xlim([l_x, u_x+5])
diagram.ax.set_ylim([l_y-squeeze_offset, u_y+squeeze_offset])

df2 = pd.DataFrame({'E': diagram.energies, 'color': diagram.colors}, diagram.positions)
for i, (k, v) in enumerate(df2.groupby(level=0)):
	v = v.sort_values('E')

	if draw_e_box:
		x = (i+1)*(diagram.dimension+diagram.space)-diagram.space/2
		y = v['E'].mean()
		w = e_box_width
		h = len(v)*e_box_offset+0.5

		if i<1: h=e_box_offset+0.5

		fancy = FancyBboxPatch((x-w/2, y-h/2), w, h, 'square', facecolor=e_box_background, edgecolor='k', alpha=0.8, zorder=3)
		diagram.ax.add_patch(fancy)
		for j, (pos, row) in enumerate(v.iterrows()):
			if i<1:
				diagram.ax.text(x, y-0.3, 0.0, color='k', ha='center', va='center', size=e_box_textsize, weight=600)
				continue

			e, c = row['E'], row['color']
			diagram.ax.text(x, y+e_box_offset*(j-(len(v)-1)/2)-0.3, round(float(e), 1), color=c, ha='center', va='center', size=e_box_textsize, weight=600)

	if draw_labels:
		x = i*(diagram.dimension+diagram.space)+diagram.dimension/2
		y = v['E'].min() - label_offset
		diagram.ax.text(x, y, name_labels[i], ha='center', size=label_textsize, weight='bold')

	if draw_images:
		x = i*(diagram.dimension+diagram.space)+diagram.dimension/2
		y = (v['E'].min() if image_pos[i] < 0 else v['E'].max())
		w, h = 15, 15

		imagebox = OffsetImage(images[i], interpolation='gaussian', zoom=image_size[i])
		ab = AnnotationBbox(imagebox, xy=(x, y+image_pos[i]), xycoords='data', frameon=False)
		ab.set_zorder(2)
		diagram.ax.add_artist(ab)

if export_png: plt.savefig(f'{output_name}.png', transparent=False, dpi=100, format="png")
if export_svg: plt.savefig(f'{output_name}.svg', transparent=True, dpi=100, format="svg")
if export_pdf: plt.savefig(f'{output_name}.pdf', transparent=True, dpi=100, format="pdf")
if show_plot: plt.show()
