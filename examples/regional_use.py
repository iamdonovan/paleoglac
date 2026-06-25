"""
Sample script to demonstrate basic use of PaleoGlac object, from initializing to calculating ELA,
paleo surface, and volume change.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import geopandas as gpd
import geoutils as gu
from paleoglac import examples, PaleoGlac


outlines = gu.Vector(examples.get_path('outlines'))
outlines.ds.set_index('lia_id', inplace=True)

elevation_changes = []
equilibrium_lines = []

for ind in outlines.index:
    glac = outlines.ds.loc[[ind]]

    # initialize the PaleoGlac and reproject to EPSG:6393
    pglac = PaleoGlac(examples.get_path('ref_dem'), glac.to_crs(6393), buffer_width=100, res=30)

    # calculate the ELA using the AABR method and a BR value of 1.56
    pglac.get_ela(method='aabr', br_val=1.56)
    outlines.ds.loc[ind, 'ela'] = pglac.ela

    # reconstruct the surface using the default ('carrivick' surface)
    pglac.reconstruct_surf()

    # append the elevation change raster and equilibrium line .ds attribute to the regional list
    elevation_changes.append(pglac.elev_change)
    equilibrium_lines.append(pglac.equilibrium_line.ds)

    # calculate the volume change in km³
    outlines.ds.loc[ind, 'volume_change'] = pglac.volume_change / 1e9

# combine the elevation changes and equilibrium lines into a single object
elevation_changes = gu.raster.merge_rasters(elevation_changes)
equilibrium_lines = gu.Vector(gpd.pd.concat(equilibrium_lines))

fig, ax = plt.subplots(1, 1, figsize=(6, 6))

outlines.to_crs(6393).boundary.plot(color='k', ax=ax, label='glacier outlines')
equilibrium_lines.plot(color='m', ax=ax, label='equilibrium lines')
elevation_changes.plot(cmap='Reds', vmax=200, vmin=0, cbar_title='elevation loss (m)')

ax.legend(loc='upper left')

ax.set_xticks([])
ax.set_yticks([])

ax.annotate(f"total volume loss: {outlines['volume_change'].sum():.2f} km$^3$",
            (0.02, 0.02),
            xycoords='axes fraction')

fig.savefig(Path('examples') / 'regional_example.png', bbox_inches='tight', dpi=200)
