import geoutils as gu
import numpy as np
from shapely.geometry import MultiPoint
from skimage.morphology import closing, disk
from skimage.filters import gaussian
from metpy import interpolate

from typing import Tuple
from geoutils.raster import RasterType
from geoutils.vector import VectorType


available_methods = [
    'carrivick_surface'
]

def ablation_area(
        dem: RasterType,
        glac: VectorType,
        ela: float | int
) -> VectorType:
    """
    Return a (simplified) representation of a glacier ablation area, defined as all parts of the glacier that are
    below the ELA:

     ablation_area = (elevation < ELA) & (onglacier)

    After identifying all DEM pixels that meet the above criteria, the mask is generalized using a morphological
    closing with a disk of radius 2 (skimage.morphology.disk), then polygonized, smoothed using a ±buffer of the DEM
    resolution, and clipped to the glacier geometry.

    :param dem: The DEM to use to determine the ablation area. Ideally the same as was used to calculate the ELA.
    :param glac: The glacier outline to use to determine on-glacier DEM pixels.
    :param ela: The glacier ELA.
    :return: the glacier ablation mask, with properties 'area' (area in km²) and 'ela' (ELA used to create the mask).
    """
    glac_buff = gu.Vector(glac).buffer(np.max(dem.res)).create_mask(dem)

    mask = (glac_buff & (dem < ela))

    closed_mask = closing(mask.data.data, footprint=disk(2))

    mask = mask.copy(new_array=(closed_mask & glac_buff).data).polygonize()
    mask = (mask.buffer(np.max(dem.res), join_style='bevel')
                .buffer(-np.max(dem.res), join_style='bevel')
                .union_all().clip(glac.union_all().ds)
            ).explode()

    mask['ela'] = ela
    mask['area'] = mask.geometry.area / 1e6

    return mask

# TODO: use/test adaptive resolution to help cut down processing time
def carrivick_surface(
        dem: RasterType,
        ablation_area: VectorType,
        simplify: bool = True
) -> RasterType:
    """
    Reconstruct the paleo glacier surface using the method described in, e.g., Carrivick et al., 2023,
    https://doi.org/10.1029/2023GL103950

    In short:

        1. Ablation area surface is interpolated using natural neighbor interpolation from elevations on
           the boundary of the ablation area surface.
        2. The difference between the DEM and the interpolated surface is calculated. Differences less than zero
           are set to 0, and the difference is smoothed using a gaussian filter with sigma=2.
        3. The smoothed difference is added to the original DEM surface.

    If simplify == True, the ablation area geometry is smoothed using a buffer of ±2.5 times the DEM pixel spacing using
    a 'bevel' join_style, before being simplified using gu.Vector.simplify with a tolerance of 0.1 times the DEM pixel
    spacing.

    :param dem: the DEM surface to use for interpolation.
    :param ablation_area: the glacier ablation area.
    :param simplify: whether to simplify the ablation area boundary.
    :return: the reconstructed ablation area ice surface.
    """
    surf_mask = ablation_area.create_mask(dem)
    tol = min(np.abs(dem.res))

    if simplify:
        boundary = (ablation_area
                    .buffer(tol*2.5, join_style='bevel')
                    .buffer(-tol*2.5, join_style='bevel')
                    .boundary.simplify(tolerance=tol/10))
    else:
        boundary = ablation_area.boundary

    pt_geoms = []
    for geom in boundary.geometry:
        if not hasattr(geom, 'geoms'):
            pt_geoms.append(MultiPoint(geom.coords))
        else:
            for geo in geom.geoms:
                pt_geoms.append(MultiPoint(geo.coords))

    xx = [pt.x for mp in pt_geoms for pt in mp.geoms]
    yy = [pt.y for mp in pt_geoms for pt in mp.geoms]

    # reduce using the maximum value within 3 pixels of each boundary point
    pc = dem.reduce_points((xx, yy), window=3, reducer_function=np.nanmax)

    XX, YY = dem.coords(force_offset='center')

    xy = np.array(list(zip(pc.geometry.x, pc.geometry.y)))
    XY = np.array(list(zip(XX[surf_mask.data], YY[surf_mask.data])))
    #ZZ = interpolate.natural_neighbor_to_grid(pc.geometry.x, pc.geometry.y, pc['z'].values, XX, YY)
    ZZ = interpolate.natural_neighbor_to_points(xy, pc['z'].values, XY)

    interp_surf = np.zeros_like(XX)
    interp_surf[surf_mask.data] = ZZ
    # interp_surf = pc.grid(dem, resampling='cubic', dist_nodata_pixel=max(dem.width, dem.height))

    interp_surf = gu.Raster(dem.copy(new_array=interp_surf))
    interp_surf.set_mask(~surf_mask)

    diff = interp_surf.data.data - dem.data.data
    diff[np.logical_or(np.isnan(diff), diff < 0)] = 0

    elev_diff = dem.copy(new_array=gaussian(diff, 2))
    elev_diff.set_mask(~surf_mask)

    return interp_surf + elev_diff
