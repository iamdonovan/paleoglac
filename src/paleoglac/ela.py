"""
Tools for calculating glacier ELA
"""
import numpy as np
import pandas as pd
from geoutils.raster import RasterType
from geoutils.vector import VectorType


np.seterr(all=None, divide=None, over=None, under=None, invalid="ignore")

def get_3d_surface_area(dem: RasterType) -> RasterType:
    """
    Approximate the 3D surface area of a DEM pixel as the product of the length of
    the slope in both the x and y direction.

    :param dem: the DEM of the surface.
    :return: a Raster of pixel area values.
    """
    dh_x = np.diff(dem.data.data, prepend=0)
    dh_y = np.diff(dem.data.data, axis=0, prepend=0)

    hyp_x = np.sqrt(dh_x**2 + dem.res[0]**2)
    hyp_y = np.sqrt(dh_y**2 + dem.res[1]**2)

    return dem.copy(new_array=(hyp_x * hyp_y))


def area_altitude_distribution(
        dem: RasterType,
        glacier_mask: VectorType,
        interval: float | int = 50,
        min_bins: int = 5,
        surface: bool = False,
        pretty_bands: bool = True,
        is_km2: bool = True
) -> pd.DataFrame:
    """
    Return a pandas DataFrame of the area altitude distribution (AAD) of a glacier or glaciers, using all DEM pixels
    that fall within the provided glacier mask.

    :param dem: the DEM to use to calculate the AAD
    :param glacier_mask: a vector representation of the glacier(s)
    :param interval: the elevation bin interval to use. Defaults to 50 m.
    :param min_bins: the minimum number of bins to calculate. If the chosen interval returns fewer than min_bins
        for the glacier's elevation range, this is recalculated to provide the minimum number of bins.
    :param surface: use an approximation of the 3D surface area, rather than the map area.
    :param pretty_bands: shift the band boundaries to be multiples of interval, rather than starting at the minimum
       glacier elevation
    :param is_km2: return the area values in units of km² [divided by 1e6]
    :return: the area altitude distribution as a DataFrame, with index equal to the bin elevation
    """

    if pretty_bands:
        min_el = max(0, interval * np.floor(dem[glacier_mask].min() / interval))
        max_el = interval * np.ceil(dem[glacier_mask].max() / interval)
    else:
        min_el = dem[glacier_mask].min()
        max_el = dem[glacier_mask].max() + interval + 1

    bins = np.arange(min_el, max_el + 1, interval)

    if len(bins) < min_bins:
        print(f"Too few bins with interval of {interval}.", end='')
        interval = np.floor((max_el - min_el) / min_bins)
        print(f" Using interval of {interval} to give {min_bins} bins.")

        bins = np.arange(min_el, max_el + 1, interval)

    if surface:
        surface_area = get_3d_surface_area(dem)
        binned_els = np.digitize(dem.data, bins)

        areas = []
        for ind, el in enumerate(bins, start=1):
            _mask = np.logical_and(binned_els == ind, glacier_mask.data)
            if np.count_nonzero(_mask) > 0:
                areas.append(surface_area[_mask].sum())
            else:
                areas.append(0)

        areas = np.array(areas)
        areas = areas[:-1]
    else:
        areas, bins = np.histogram(dem[glacier_mask], bins)

        areas = areas.astype(float)
        areas *= dem.res[0] * dem.res[1]

    if is_km2:
        areas /= 1e6

    return pd.DataFrame(data={'elevation': bins[:-1],
                              'elevation_mid': bins[:-1] + interval / 2,
                              'area': areas}).set_index('elevation')


def mge_ela(
        dem: RasterType,
        glac_mask: RasterType,
        interval: float | int = 50
) -> float:
    """
    Calculate the glacier ELA using the Median Glacier Elevation (MGE) method (Kurowski, 1891; Sissons, 1974).

    Here, the ELA is calculated as the area-weighted median elevation of the glacier.

    :param dem: the DEM of the glacier surface
    :param glac_mask: the glacier mask
    :param interval: the contour interval to use to calculate the AAD.
    :return: the MGE ELA value.
    """
    el_df = area_altitude_distribution(dem, glac_mask, interval, pretty_bands=False, surface=True)
    aar_vals = el_df['area'] * el_df['elevation_mid']

    return aar_vals.sum() / el_df['area'].sum()


def aabr_ela(
        dem: RasterType,
        glac_mask: RasterType,
        interval: float | int = 50,
        br_val: float = 1.75
) -> float:
    """
    Calculate the ELA using the Area-Altitude Balance Ratio (AABR) method
    (e.g., Osmaston, 2005, https://doi.org/10.1016/j.quaint.2005.02.004).

    :param dem: the DEM of the glacier surface
    :param glac_mask: the glacier mask
    :param interval: the contour interval to use to calculate the AAD.
    :param br_val: the balance ratio (BR) value to use. Default is 1.75, the "global" value suggested by
        Rea (2009), https://doi.org/10.1016/j.quascirev.2008.10.011.
    :return: the AABR ELA value
    """
    el_df = area_altitude_distribution(dem, glac_mask, interval, pretty_bands=False, surface=True).reset_index()

    for ind, band in el_df.iterrows():
        aabr_vals = el_df['area'] * (el_df['elevation_mid'] - band.elevation)
        aabr_vals[aabr_vals < 0] *= br_val

        el_df.loc[ind, 'aabr_sum'] = aabr_vals.sum()

    low_ind = el_df.loc[el_df['aabr_sum'] < 0].index[0]

    interp_ela = el_df.loc[low_ind, 'elevation'] - interval
    interp_ela += (interval * abs(el_df.loc[low_ind - 1, 'aabr_sum'])) / (
                abs(el_df.loc[low_ind - 1, 'aabr_sum']) + abs(el_df.loc[low_ind, 'aabr_sum']))

    return int(interp_ela)


def aar_ela(
        dem: RasterType,
        glac_mask: RasterType,
        interval: float | int = 50,
        aar_val: float = 0.56
) -> float:
    """
    Calculate the ELA as the elevation whose cumulative area, starting from the lowest elevation,
        is closest to (1 - aar_val) * total_area.

    :param dem: the DEM of the glacier surface
    :param glac_mask: the glacier mask
    :param interval: the contour interval to use to calculate the AAD.
    :param aar_val: the proposed accumulation area ratio (AAR). Default is 0.56, based on WGI-derived values from
        Kern and László (2010), https://doi.org/10.1016/j.quascirev.2010.06.033
    :return: the AAR ELA
    """
    aad = area_altitude_distribution(dem, glac_mask, interval, pretty_bands=False, surface=True)
    accumulation_area = aar_val * aad.area.sum()

    aad['cumulative_area'] = aad.area.cumsum()

    # get the elevations of the bands closest to the accumulation area boundary
    low = aad.loc[aad['cumulative_area'] < accumulation_area, 'elevation_mid'].index[-1]
    high = aad.loc[aad['cumulative_area'] > accumulation_area, 'elevation_mid'].index[0]

    # get the proportion of area below the accumulation area in these two bands
    frac = (accumulation_area - aad.loc[low, 'cumulative_area']) / (
                aad.loc[high, 'cumulative_area'] - aad.loc[low, 'cumulative_area'])

    return low + frac * interval
