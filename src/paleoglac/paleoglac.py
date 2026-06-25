import pathlib
from collections import abc

import pandas as pd
import geopandas as gpd

import rasterio as rio
from rasterio.crs import CRS
from rasterio.enums import Resampling

import geoutils as gu
from geoutils._typing import DTypeLike
from geoutils.raster import RasterType
from geoutils.raster.distributed_computing.multiproc import MultiprocConfig
from geoutils.vector import VectorType

from typing import Any, Literal

from . import ela, surface


glac_attrs = ('ela',
              'equilibrium_line'
              'ablation_area',
              'ablation_area_mask',
              'surface',
              'bed',
              'elev_change')

class PaleoGlac(gu.Raster):

    def __init__(
            self,
            filename_or_dataset: (
                str | pathlib.Path | RasterType | rio.io.DatasetReader | rio.io.MemoryFile | dict[str, Any] | PaleoGlac
            ),
            geom: VectorType | gpd.GeoSeries | gpd.GeoDataFrame = None,
            buffer_width: int = 0,
            crs: CRS | str | int | None = None,
            res: float | abc.Iterable[float] | None = None,
            crop_to_geom: bool = True
    ) -> None:
        """
        A representation of a paleoglacier (or, alternatively, an extant glacier), with the following basic attributes:

            - data: the surface elevation, a RasterType
            - glacier_outline: the bounds of the glacier, a VectorType
            - glacier_mask: a binary mask of glacier/not glacier pixels, based on the outline

        Additional attributes can be set/created after initialization:

            - ela: the glacier equilibrium line altitude. For details on implemented/available methods,
              see PaleoGlac.get_ela
            - ablation_area: a vector representation of the on-glacier area below the calculated ELA
            - ablation_area_mask: a raster mask of the ablation area
            - equilibrium_line: a vector representation of the equilibrium line
            - aar: the glacier accumulation area ratio
            - paleo_surface: the reconstructed paleo surface. For details on implemented methods, see
              PaleoGlac.reconstruct_surface
            - elev_change: the elevation difference between the current and reconstructed surface
            - volume_change: the volume difference between the current and reconstructed surface
            - mean_thick_change: the mean thickness change (i.e., mean of all elevation changes)

        :param filename_or_dataset: Path to file or Rasterio dataset representing the glacier surface elevation
        :param geom: the glacier geometry, in the form of a geoutils Vector or geopandas GeoSeries / GeoDataFrame
        :param buffer_width: the distance (in the units of the CRS) to buffer the geometry by before cropping the DEM
        :param crs: Coordinate reference system. Any CRS supported by Pyproj (e.g., CRS object, EPSG integer).
        :param res: Destination resolution (pixel size) in units of destination CRS. Single value or (xres, yres).
        :param crop_to_geom: crop the input DEM to the glacier geometry
        """
        super().__init__(filename_or_dataset)
        dem_crs = self.crs

        # figure out what crs we should use
        # priority: crs -> geom.crs -> dem.crs
        if crs is not None:
            out_crs = crs
        elif geom is not None:
            out_crs = geom.crs
        else:
            out_crs = dem_crs

        # initialize empty values if they don't already exist
        if not isinstance(filename_or_dataset, PaleoGlac):
            self._ela = None
            self._ablation_area = None
            self._ablation_area_mask = None
            self._equilibrium_line = None
            self._aar = None
            self._paleo_surface = None

        if geom is not None:
            self.glacier_outline = gu.Vector(geom)

            # only crop if asked to crop
            if crop_to_geom:
                self.crop(self.glacier_outline.to_crs(out_crs).buffer(buffer_width), inplace=True)

            # only reproject if there's a CRS mismatch, or it's been requested
            if self.crs != out_crs or res is not None:
                self.reproject(crs=out_crs, res=res, inplace=True)

            # add the area of the outline in the current crs

        else:
            # re-project the outline and any other geometries if needed?
            if self.crs != out_crs and res is not None:
                self.reproject(crs=out_crs, res=res, inplace=True)


    @property
    def glacier_outline(self) -> gu.Vector:
        return self._glacier_outline

    @glacier_outline.setter
    def glacier_outline(
            self: PaleoGlac,
            geom: gu.Vector | gpd.GeoSeries | gpd.GeoDataFrame) -> None:
        if isinstance(geom, gu.Vector):
            self._glacier_outline = geom
        elif isinstance(geom, (gpd.GeoSeries, gpd.GeoDataFrame)):
            self._glacier_outline = gu.Vector(geom)
        else:
            raise ValueError("The glacier geometry must be a Vector, GeoSeries, or GeoDataFrame.")

    @property
    def glacier_mask(self: PaleoGlac):
        return self._glacier_outline.create_mask(self)

    @property
    def boundary(self: PaleoGlac) -> VectorType:
        return self.glacier_outline.boundary

    @property
    def ela(self: PaleoGlac) -> float:
        """
        The calculated ELA for the glacier.
        """
        return self._ela

    @ela.setter
    def ela(self: PaleoGlac,
            value: float | int) -> None:
        self._ela = value

    @property
    def ablation_area(self: PaleoGlac) -> VectorType:
        if self._ela is None:
            raise ValueError("ELA must be set using .get_ela() before calculating ablation area.")
        return surface.ablation_area(self, self.glacier_outline, self.ela)

    @property
    def ablation_area_mask(self: PaleoGlac):
        return self.ablation_area.create_mask(self)

    @property
    def equilibrium_line(self: PaleoGlac) -> VectorType:
        """
        The glacier's equilibrium line, calculated as the difference between the glacier
        boundary and ablation area boundary.
        """
        return self.ablation_area.boundary.overlay(self.boundary, how='difference')

    @property
    def aar(self: PaleoGlac) -> float:
        """
        the glacier Accumulation Area Ratio
        """
        total_area = self._glacier_outline.area.sum()
        abl_area = self.ablation_area.area.sum()
        return (total_area - abl_area) / total_area

    @property
    def paleo_surface(self: PaleoGlac) -> RasterType:
        return self._paleo_surface

    @paleo_surface.setter
    def paleo_surface(self: PaleoGlac, surf: RasterType) -> None:
        self._paleo_surface = surf

    @property
    def elev_change(self: PaleoGlac) -> RasterType:
        if self._paleo_surface is not None:
            return self._paleo_surface - self
        else:
            raise ValueError("Surface must be reconstructed before calculating elevation change.")

    @property
    def volume_change(self: PaleoGlac) -> RasterType:
        return self.elev_change.data.sum() * self.res[0] * self.res[1]

    @property
    def mean_thick_change(self: PaleoGlac) -> float:
        return self.elev_change.data.mean()

    def reproject(
            self: PaleoGlac,
            ref: RasterType | str | None = None,
            crs: CRS | str | int | None = None,
            res: float | abc.Iterable[float] | None = None,
            grid_size: tuple[int, int] | None = None,
            bounds: dict[str, float] | rio.coords.BoundingBox | None = None,
            nodata: int | float | None = None,
            dtype: DTypeLike | None = None,
            resampling: Resampling | str = Resampling.bilinear,
            force_source_nodata: int | float | None = None,
            inplace: bool = False,
            silent: bool = False,
            n_threads: int = 0,
            memory_limit: int = 64,
            multiproc_config: MultiprocConfig | None = None,
    ) -> RasterType | None:

        inargs = locals()
        inargs.pop('self')
        inargs.pop('__class__')

        # figure out what crs we should use
        # priority: crs -> ref.crs -> dem.crs
        if crs is not None:
            pass
        elif ref is not None:
            crs = ref.crs
        else:
            crs = self.crs

        if not inplace:
            out_obj = super().reproject(**inargs)

            out_obj.glacier_outline = self.glacier_outline
            out_obj.update_outline_masks(crs)

            # if the ela is already set, copy it over
            out_obj.ela = self.ela

            if self.paleo_surface is not None:
                out_obj.paleo_surface = self.paleo_surface.reproject(out_obj)

            return out_obj

        else:
            super().reproject(**inargs)
            self.update_outline_masks(crs)

            if self.paleo_surface is not None:
                self.paleo_surface = self.paleo_surface.reproject(self)

    def crop(
            self: PaleoGlac,
            bbox: PaleoGlac | RasterType | gu.Vector | list[float] | tuple[float, ...],
            mode: Literal["match_pixel"] | Literal["match_extent"] = "match_pixel",
            *,
            inplace: Literal[False] = False
    ) -> PaleoGlac | RasterType | None:

        inargs = locals()
        inargs.pop('self')
        inargs.pop('__class__')

        if not inplace:
            out_obj = super().crop(**inargs)
            out_obj.glacier_outline = self.glacier_outline

            out_obj.update_outline_masks(None)

            # if the ela is already set, copy it over
            out_obj.ela = self.ela

            if self.paleo_surface is not None:
                out_obj.paleo_surface = self.paleo_surface.crop(out_obj)

            return out_obj
        else:
            super().crop(**inargs)
            self.update_outline_masks(None)
            if self.paleo_surface is not None:
                self.paleo_surface = self.paleo_surface.crop(self)

            return None

    def update_outline_masks(
            self: PaleoGlac,
            crs: CRS | str | int | None = None) -> None:
        """
        Update / reproject outlines and masks after reprojecting.

        :param crs: Destination coordinate reference system as a string or EPSG. If ``ref`` not set,
            defaults to this raster's CRS.
        """
        if crs is not None:
            self.glacier_outline = self.glacier_outline.to_crs(crs)

    # different ways of interpolating the surface
    def reconstruct_surf(
            self: PaleoGlac,
            set_val: bool = True,
            method: Literal['carrivick'] ='carrivick',
            **kwargs) -> RasterType | None:
        """
        Reconstruct the paleo surface using the given methods. For a list of available methods, see
        paleoglac.surface.available_methods.

        :param set_val: set the interpolated surface as object's .paleo_surface property.
        :param method: the name of the method to use for interpolating the paleo surface.
        :param kwargs: additional kwargs for the chosen method.
        :return: the paleo surface, if set_val is False
        """
        #method determines which algorithm we use
        #methods: carrivick, ...
        methods = {
            'carrivick': surface.carrivick_surface,
        }

        surf = methods[method](self, self.ablation_area, **kwargs)

        if set_val:
            self.paleo_surface = surf
            return None
        else:
            return surf

    # implement different methods of getting the ELA
    def get_ela(
            self: PaleoGlac,
            set_val: bool = True,
            method: Literal['aabr', 'aar', 'mge', 'kurowski'] = 'aabr',
            interval: int | float = 50,
            **kwargs) -> float:
        """
        Calculate the glacier ELA using the requested method. Available methods are:

            - aabr (e.g., Osmaston, 2005, https://doi.org/10.1016/j.quaint.2005.02.004) [default]
            - aar (e.g., Porter, 2000, https://doi.org/10.1016/S0277-3791(00)00178-5)
            - mge (e.g., Sissons, 1974, https://doi.org/10.2307/621517)

        For additional details, including default parameters, see paleoglac.ela.{method}_ela.

        :param set_val: set the resulting ELA value as this object's .ela value.
        :param method: the name of the method to use for calculating the ELA
        :param interval: the elevation bin interval to use. Defaults to 50 m.
        :param kwargs: additional kwargs for the chosen method.
        :return:
        """
        #method determines which algorithm we use
        #methods: aar, aabr, mge/kurowski
        methods = {
            'aar': ela.aar_ela,
            'aabr': ela.aabr_ela,
            'mge': ela.mge_ela
        }

        ela_val = methods[method](
            self,
            self.glacier_mask,
            interval=interval,
            **kwargs
        )

        if set_val:
            self.ela = ela_val
        else:
            return ela_val

        # calculate using the surface dem
    def area_altitude_distribution(self,
                      interval: float | int = 50,
                      min_bins: int = 5,
                      surface: bool = False,
                      pretty_bands: bool = False,
                      is_km2: bool = True) -> pd.DataFrame:
        """
        Return a pandas DataFrame of the area altitude distribution of the glacier, using all DEM pixels that fall
        within the glacier's outline.

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
        return ela.area_altitude_distribution(
            self,
            self.glacier_mask,
            interval=interval,
            min_bins=min_bins,
            surface=surface,
            pretty_bands=pretty_bands,
            is_km2=is_km2
        )
