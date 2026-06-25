"""
Sample script to demonstrate basic use of PaleoGlac object, from initializing to calculating ELA,
paleo surface, and volume change.
"""
from pathlib import Path
import geoutils as gu
from paleoglac import examples, PaleoGlac
import matplotlib.pyplot as plt


outlines = gu.Vector(examples.get_path('outlines'))

# choose a single glacier to use
ind = 'LIA-01.2134'
glac = outlines.ds.loc[outlines['lia_id'] == ind]

# initialize the glacier, reprojected to EPSG 6393
pglac = PaleoGlac(examples.get_path('ref_dem'), glac.to_crs(6393), buffer_width=100, res=30)

# calculate ela using AABR method, BR of 1.56, setting as .ela property
pglac.get_ela(method='aabr', br_val=1.56)

# plot boundary and equilibrium line
fig, ax = plt.subplots(1, 1)

pglac.plot(add_cbar=False, cmap='gray', ax=ax)
pglac.boundary.plot(color='b', ax=ax, label='glacier outline')
pglac.equilibrium_line.plot(color='m', ax=ax, label=f"eq. line (ELA={pglac.ela} m)")

ax.legend(loc='upper left')

ax.set_xticks([])
ax.set_yticks([])

# reconstruct the surface using the default ('carrivick' surface)
pglac.reconstruct_surf()

# plot the elevation change
pglac.elev_change.plot(vmax=200, vmin=0, cmap='Reds', cbar_title='elevation loss (m)')

# annotate the volume change (in km³)
ax.annotate(f"total volume loss: {pglac.volume_change / 1e9:.2f} km$^3$", (0.02, 0.02), xycoords='axes fraction')

fig.savefig(Path('examples') / 'basic_example.png', bbox_inches='tight', dpi=200)

