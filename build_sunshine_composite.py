"""
The script does the following:

- Crop CM SAF SARAH-3 SDU (Sunshine Duration) NetCDF files to a bounding box
- Average multiple monthly files into a representative annual sunshine grid
- Write it as a self-describing GeoTIFF composite
- Generate a JSON colormap sidecar (dynamic grey-to-orange range, 
  target reprojection info, real image bounds, and 
  the client-side zoom-resampling hint for MapLibre GL)
- Render a PNG preview using exactly that sidecar

Each source NetCDF's scale_factor/add_offset/units are read and validated
per file before any pixel value is used: SDU values are kept packed as
Int16 raw digital numbers and must be scaled by the provided factor to
recover physical hours (physical hours = raw * scale_factor + add_offset).

The sidecar's `image_bounds` field is the single source of truth for where
the PNG belongs on a map: it is read directly from the same edges used to
georeference the GeoTIFF, used by client code (e.g. a MapLibre ImageSource).

Cropping: pixel windows and their geographic edges are derived directly
from the SDU band's own affine transform via
`rasterio.windows.from_bounds`/`rasterio.windows.bounds`. GDAL guarantees
a dataset's `.read()` output and its own `.transform` are always mutually
consistent -- that is the definition of a valid raster dataset -- so
anchoring on one dataset's own transform makes a coordinate/content
mismatch structurally impossible.

Reprojection note: the PNG is rendered into Web Mercator (EPSG:3857) by
default. This composite spans roughly 30-72 deg N -- wide enough that a
map image overlay needs its source pixels already laid out linearly in
Mercator-y for the raster to align correctly with a Mercator-based
basemap across its full extent, not just at its four corners.
"""

import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image
from rasterio.crs import CRS
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, calculate_default_transform, reproject
from rasterio.windows import Window, bounds as window_bounds, from_bounds as window_from_bounds

SOURCE_CRS = CRS.from_epsg(4326)
WEBMERCATOR_CRS = CRS.from_epsg(3857)
FILL_VALUE = -1.0
NODATA_VALUE = -9999.0
EXPECTED_UNIT = "h"

# Matches: SDUms YYYYMMDD HHMMSS ... e.g. SDUms2021010100000042310001I1MA.nc
FILENAME_DATE_RE = re.compile(r"^SDUms(\d{4})(\d{2})\d{2}\d{6}")


@dataclass(frozen=True)
class BoundingBox:
    west: float
    south: float
    east: float
    north: float


@dataclass(frozen=True)
class YearMonth:
    year: int
    month: int

    def __le__(self, other: "YearMonth") -> bool:
        return (self.year, self.month) <= (other.year, other.month)

    def __ge__(self, other: "YearMonth") -> bool:
        return (self.year, self.month) >= (other.year, other.month)


@dataclass(frozen=True)
class ScaleOffset:
    scale: float
    offset: float
    unit: str
    nodata: float


def parse_year_month(path: Path) -> YearMonth:
    """Extract (year, month) from a CM SAF SDU monthly filename."""

    match = FILENAME_DATE_RE.match(path.name)

    if not match:
        raise ValueError(f"Filename does not match expected SDU pattern: {path.name}")

    year, month = match.groups()
    return YearMonth(int(year), int(month))


def discover_monthly_files(directory: Path, start: YearMonth, end: YearMonth) -> list[Path]:
    """Find SDU monthly NetCDF files in directory within [start, end] inclusive, sorted chronologically."""

    candidates = sorted(directory.glob("SDUms*.nc"))
    selected = [f for f in candidates if start <= parse_year_month(f) <= end]

    if not selected:
        raise ValueError(f"No SDU files found in {directory} between {start} and {end}.")

    return selected


def resolve_sdu_uri(nc_path: Path) -> str:
    """Find the SDU subdataset URI inside a CM SAF NetCDF file."""

    with rasterio.open(nc_path) as parent:
        subs = parent.subdatasets

    return next(s for s in subs if s.endswith(":SDU"))


def read_scale_offset(sdu_uri: str) -> ScaleOffset:
    """
    Read the SDU band's scale_factor, add_offset, units, and nodata directly from the file.

    These are read per-file rather than assumed, since a corrupted download,
    a different product version, or a mixed-in unrelated file could carry
    different packing parameters.
    """

    with rasterio.open(sdu_uri) as src:

        scale = src.scales[0] if src.scales else 1.0
        offset = src.offsets[0] if src.offsets else 0.0
        unit = (src.units[0] if src.units else None) or src.tags(1).get("units", "")

        nodata = src.nodata

        if nodata is None:
            nodata = float(src.tags(1).get("_FillValue", "nan"))

    return ScaleOffset(scale=float(scale), offset=float(offset), unit=str(unit), nodata=float(nodata))


def validate_scale_offset(so: ScaleOffset, source_path: Path) -> None:
    """
    Raise ValueError if the scale/offset/unit read from a file are missing or physically impossible.

    A scale_factor must be a finite, strictly positive number: zero would
    collapse all data to the offset, a negative value would invert the
    physical quantity, and NaN/inf indicates a missing or corrupted attribute.
    """

    if so.scale is None or not math.isfinite(so.scale):
        raise ValueError(f"Invalid scale_factor in {source_path.name}: {so.scale!r} is not a finite number.")

    if so.scale <= 0:
        raise ValueError(
            f"Impossible scale_factor in {source_path.name}: {so.scale} must be strictly positive "
            "for a physical duration quantity."
        )

    if so.offset is None or not math.isfinite(so.offset):
        raise ValueError(f"Invalid add_offset in {source_path.name}: {so.offset!r} is not a finite number.")

    if so.unit != EXPECTED_UNIT:
        raise ValueError(
            f"Unexpected units in {source_path.name}: got '{so.unit}', expected '{EXPECTED_UNIT}' (hours). "
            "Refusing to proceed since downstream code assumes hours."
        )

    if not math.isfinite(so.nodata):
        raise ValueError(f"Invalid nodata/_FillValue in {source_path.name}: {so.nodata!r} is not a finite number.")


def validate_consistent_scaling(files: list[Path]) -> ScaleOffset:
    """
    Read and validate scale/offset/unit for every file, and require they all agree.

    Averaging monthly grids that were packed with different scale_factor
    values would silently mix incompatible units. Returns the shared
    ScaleOffset once consistency is confirmed.
    """

    reference = None
    reference_path = None

    for path in files:
        sdu_uri = resolve_sdu_uri(path)
        so = read_scale_offset(sdu_uri)
        validate_scale_offset(so, path)

        if reference is None:
            reference, reference_path = so, path
            continue

        if not math.isclose(so.scale, reference.scale, rel_tol=1e-6):
            raise ValueError(
                f"Inconsistent scale_factor across inputs: {path.name} has {so.scale}, "
                f"but {reference_path.name} has {reference.scale}."
            )
        if not math.isclose(so.offset, reference.offset, abs_tol=1e-9):
            raise ValueError(
                f"Inconsistent add_offset across inputs: {path.name} has {so.offset}, "
                f"but {reference_path.name} has {reference.offset}."
            )
        if so.unit != reference.unit:
            raise ValueError(
                f"Inconsistent units across inputs: {path.name} has '{so.unit}', "
                f"but {reference_path.name} has '{reference.unit}'."
            )

    print(f"Validated scale_factor={reference.scale}, add_offset={reference.offset}, "
          f"units='{reference.unit}' across {len(files)} files")
    return reference


def compute_crop_window(sdu_uri: str, bbox: BoundingBox) -> Window:
    """
    Compute the pixel window covering bbox, directly from the SDU dataset's
    own affine transform (as GDAL derives it from the file's CF lon/lat
    coordinate variables).

    This deliberately avoids cross-referencing a separately-opened
    `lat_bnds`/`lon_bnds` variable to find crop indices: GDAL guarantees
    that `.read()` and `.transform` on the SAME dataset object always agree
    on pixel layout, but makes no such guarantee across two different
    variables read via separate `rasterio.open()` calls. Anchoring on one
    dataset's own transform makes a content/coordinate mismatch (like the
    real-world content landing several degrees away from where its own
    written coordinates say it is) structurally impossible.
    """
    with rasterio.open(sdu_uri) as src:
        window = window_from_bounds(
            bbox.west, bbox.south, bbox.east, bbox.north, transform=src.transform
        )
        window = window.round_offsets().round_lengths()
        full_height, full_width = src.height, src.width

    row_start = max(0, window.row_off)
    col_start = max(0, window.col_off)
    row_stop = min(full_height, window.row_off + window.height)
    col_stop = min(full_width, window.col_off + window.width)

    if row_stop <= row_start or col_stop <= col_start:
        raise ValueError("Bounding box does not intersect any pixels in the grid.")

    return Window(col_off=col_start, row_off=row_start, width=col_stop - col_start, height=row_stop - row_start)


def window_edges(sdu_uri: str, window: Window) -> BoundingBox:
    """
    Compute the true pixel-edge bounding box for a window, from the SDU
    dataset's own transform.

    This is the single authoritative source for the composite's real
    geographic extent: it is used both to build the GeoTIFF's affine
    transform and to populate the sidecar's `image_bounds`, so the two can
    never drift apart, and it is derived from the exact same transform used
    to actually read the pixel data, so content and coordinates can never
    silently disagree either.
    """
    with rasterio.open(sdu_uri) as src:
        west, south, east, north = window_bounds(window, src.transform)

    return BoundingBox(west=float(west), south=float(south), east=float(east), north=float(north))


def read_cropped_sdu(sdu_uri: str, window: Window, scale_offset: ScaleOffset) -> tuple[np.ndarray, np.ndarray]:
    """
    Read the SDU band's pixel window and unpack raw digital numbers into physical hours.

    Reading via `rasterio`'s windowed read on the SDU dataset directly
    (rather than reading the full array and then manually slicing/flipping
    it based on a separately-read coordinate variable) guarantees the
    returned array is already in the correct, north-up orientation implied
    by the dataset's own `.transform` -- no manual axis-order correction is
    needed for either dimension.

    Nodata detection happens on the raw digital number (matching the
    file's stored _FillValue) before the scale/offset conversion is
    applied, since packed nodata sentinels are not meant to be scaled.

    Returns (sdu_crop_hours, valid_mask).
    """

    with rasterio.open(sdu_uri) as src:
        raw_crop = src.read(1, window=window).astype(np.float64)

    valid = np.isfinite(raw_crop) & ~np.isclose(raw_crop, scale_offset.nodata)

    sdu_crop = (raw_crop * scale_offset.scale + scale_offset.offset).astype(np.float32)

    return sdu_crop, valid


def compute_average_annual_sdu(
    files: list[Path], window: Window, scale_offset: ScaleOffset
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build a representative average annual sunshine grid (in hours) from monthly files.

    For each calendar month (Jan..Dec), average that month's grid across all
    available years in `files`, then sum the 12 monthly averages into one
    average annual total. Returns (annual_grid, valid_mask).
    """

    by_month: dict[int, list[np.ndarray]] = defaultdict(list)
    valid_by_month: dict[int, list[np.ndarray]] = defaultdict(list)

    for path in files:

        year_month = parse_year_month(path)
        sdu_uri = resolve_sdu_uri(path)
        sdu_crop, valid = read_cropped_sdu(sdu_uri, window, scale_offset)
        by_month[year_month.month].append(sdu_crop)
        valid_by_month[year_month.month].append(valid)

    monthly_averages = []
    combined_valid = None

    for month in range(1, 13):
        if month not in by_month:
            raise ValueError(f"No data found for calendar month {month:02d} in the given file range.")

        stacked = np.stack(by_month[month])
        stacked_valid = np.stack(valid_by_month[month])

        month_valid = stacked_valid.any(axis=0)
        safe_stacked = np.where(stacked_valid, stacked, 0.0)
        month_sum = safe_stacked.sum(axis=0)
        month_count = stacked_valid.sum(axis=0)

        month_average = np.divide(
            month_sum, month_count, out=np.zeros_like(month_sum), where=month_count > 0
        )

        monthly_averages.append(month_average)
        combined_valid = month_valid if combined_valid is None else (combined_valid & month_valid)

    annual_grid = np.sum(monthly_averages, axis=0).astype(np.float32)
    return annual_grid, combined_valid


def write_annual_grid_geotiff(
    annual_grid: np.ndarray,
    valid: np.ndarray,
    edges: BoundingBox,
    files: list[Path],
    start: YearMonth,
    end: YearMonth,
    scale_offset: ScaleOffset,
    output_path: Path,
) -> tuple[float, float]:
    """
    Write the average annual SDU composite (already unpacked to hours) as a georeferenced GeoTIFF in WGS84.

    Includes: CRS + affine transform + nodata + dtype (minimal requirements),
    plus band description/units (strongly recommended additions), the
    validated source scale_factor/add_offset for provenance, and a set of
    categorical/provenance tags useful for pipeline maintenance.

    Returns (data_min, data_max) over valid pixels, for reuse by the sidecar.
    """

    transform = from_bounds(edges.west, edges.south, edges.east, edges.north,
                             annual_grid.shape[1], annual_grid.shape[0])

    data_min = float(annual_grid[valid].min())
    data_max = float(annual_grid[valid].max())

    data = np.where(valid, annual_grid, NODATA_VALUE).astype(np.float32)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(
        output_path, "w",
        driver="GTiff",
        height=data.shape[0], width=data.shape[1],
        count=1, dtype="float32",
        crs=SOURCE_CRS, transform=transform,
        nodata=NODATA_VALUE,
        compress="DEFLATE",
        tiled=True,
    ) as dst:

        dst.write(data, 1)
        dst.set_band_description(1, "Average annual sunshine duration (hours)")
        dst.update_tags(
            1,
            units=scale_offset.unit,
            RepresentationType="ATHEMATIC",
            valid_min=data_min,
            valid_max=data_max,
            source_scale_factor=scale_offset.scale,
            source_add_offset=scale_offset.offset,
            AREA_OR_POINT="Area",
            variable_id="SDU",
            source="CM SAF SARAH-3 (SDU, Heliosat Edition 3)",
            institution="EUMETSAT/CMSAF",
            time_coverage_start=f"{start.year}-{start.month:02d}-01",
            time_coverage_end=f"{end.year}-{end.month:02d}-01",
            source_file_count=str(len(files)),
            source_files=";".join(f.name for f in files),
            processing="Raw Int16 unpacked via validated scale_factor/add_offset, "
                       "monthly SDU averaged per calendar month across years, then summed to annual total",
            crop_bbox=f"west={edges.west},south={edges.south},east={edges.east},north={edges.north}",
            created=datetime.now(timezone.utc).isoformat(),
            pipeline="build_sunshine_composite.py",
        )

    print(f"Wrote GeoTIFF composite: {output_path} (min={data_min:.2f}, max={data_max:.2f}), "
          f"real bounds west={edges.west}, south={edges.south}, east={edges.east}, north={edges.north}")

    return data_min, data_max


def grey_to_orange_stops() -> list[tuple[float, tuple[int, int, int]]]:
    """
    Fixed-shape colormap: dark grey (low) -> pale -> orange (high).

    Fractions are relative positions along the *dynamic* [data_min, data_max]
    range computed per dataset; only the RGB anchor colors are fixed.
    """

    return [
        (0.00, (60, 60, 60)),      # dark charcoal grey
        (0.25, (140, 140, 140)),   # medium grey
        (0.50, (230, 220, 190)),   # pale cream / warm beige
        (0.75, (255, 170, 60)),    # light amber orange
        (1.00, (230, 81, 0)),      # deep/burnt orange (Material Design "Deep Orange 900")
    ]


def write_colormap_sidecar(
    tif_path: Path,
    data_min: float,
    data_max: float,
    nodata: float,
    edges: BoundingBox,
    target_crs: CRS,
    resampling: str,
    zoom_resampling: str,
    sidecar_path: Path,
) -> None:
    """
    Write a JSON sidecar describing the dynamic colormap, target reprojection,
    the real image bounds, and the client-side zoom-resampling hint.

    `image_bounds` is copied directly from the same `edges` used to build the
    GeoTIFF's affine transform, in plain WGS84 lon/lat -- this is exactly what
    a MapLibre GL `ImageSource.coordinates` needs, so client code should read
    it from here rather than hardcoding corner values that could silently
    drift out of sync with the pipeline's actual bbox/crop.

    `reprojection.target_crs`/`resampling` describe the warp
    render_png_from_sidecar() applies before coloring the raster. The default
    is Web Mercator, since this composite's latitude extent is wide enough
    that a MapLibre `ImageSource`'s corner-only bilinear warp needs
    Mercator-linear source pixels to line up correctly across the raster's
    interior, not just at its four corners.
    `rendering.zoom_resampling` maps onto MapLibre's `raster-resampling` paint
    property, controlling GPU filtering above native zoom ("linear" = smooth,
    "nearest" = blocky).

    None of this is canonical GeoTIFF/NetCDF metadata -- it is an explicit
    convention read by the renderer and by client code.
    """

    sidecar = {
        "source_tif": tif_path.name,
        "created": datetime.now(timezone.utc).isoformat(),
        "data_range": {"min": data_min, "max": data_max},
        "nodata": nodata,
        "image_bounds": {
            "west": edges.west,
            "south": edges.south,
            "east": edges.east,
            "north": edges.north,
        },
        "colormap": {
            "name": "grey_to_orange_dynamic",
            "space": "RGB",
            "stops": [{"fraction": frac, "rgb": list(rgb)} for frac, rgb in grey_to_orange_stops()],
        },
        "reprojection": {
            "target_crs": f"EPSG:{target_crs.to_epsg()}",
            "resampling": resampling,
        },
        "rendering": {
            "zoom_resampling": zoom_resampling,
        },
    }

    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    with open(sidecar_path, "w") as f:
        json.dump(sidecar, f, indent=2)

    print(f"Wrote colormap sidecar: {sidecar_path}")


def apply_colormap(values: np.ndarray, valid: np.ndarray, data_min: float, data_max: float,
                    stops: list[tuple[float, tuple[int, int, int]]]) -> np.ndarray:
    """Map normalized values to RGBA using piecewise-linear interpolation between stops."""

    fractions = np.zeros_like(values, dtype=np.float32)
    if data_max > data_min:
        fractions[valid] = (values[valid] - data_min) / (data_max - data_min)

    fractions = np.clip(fractions, 0.0, 1.0)

    stop_fracs = np.array([s[0] for s in stops], dtype=np.float32)
    stop_rgb = np.array([s[1] for s in stops], dtype=np.float32)

    rgba = np.zeros((*values.shape, 4), dtype=np.uint8)

    for channel in range(3):
        rgba[..., channel] = np.clip(
            np.interp(fractions, stop_fracs, stop_rgb[:, channel]), 0, 255
        ).astype(np.uint8)
    rgba[..., 3] = np.where(valid, 255, 0).astype(np.uint8)

    return rgba


def render_png_from_sidecar(tif_path: Path, sidecar_path: Path, output_png_path: Path) -> None:
    """
    Render a PNG preview strictly from the GeoTIFF + its colormap sidecar.

    Reprojects to the sidecar's target CRS (Web Mercator by default), applies
    the sidecar's dynamic grey-to-orange colormap over the sidecar's stored
    data range, and saves. The reprojection here does NOT change the
    geographic footprint -- an axis-aligned WGS84 rectangle reprojects to an
    axis-aligned Web Mercator rectangle, so `image_bounds` in the sidecar
    (always stored in WGS84) remains valid for placing this PNG on a map via
    a MapLibre `ImageSource`, regardless of the PNG's own pixel CRS.
    """

    with open(sidecar_path) as f:
        sidecar = json.load(f)

    data_min = sidecar["data_range"]["min"]
    data_max = sidecar["data_range"]["max"]
    nodata = sidecar["nodata"]

    target_crs = CRS.from_string(sidecar["reprojection"]["target_crs"])
    resampling_name = sidecar["reprojection"]["resampling"]

    stops = [(s["fraction"], tuple(s["rgb"])) for s in sidecar["colormap"]["stops"]]

    with rasterio.open(tif_path) as src:

        src_data = src.read(1)
        src_crs = src.crs
        src_transform = src.transform
        src_nodata = src.nodata if src.nodata is not None else nodata

        if target_crs == src_crs:
            destination = src_data
            valid_dst = destination != src_nodata
        else:
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src_crs, target_crs, src_data.shape[1], src_data.shape[0],
                left=src_transform.c, top=src_transform.f,
                right=src_transform.c + src_transform.a * src_data.shape[1],
                bottom=src_transform.f + src_transform.e * src_data.shape[0],
            )

            destination = np.full((dst_height, dst_width), src_nodata, dtype=np.float32)
            reproject(
                source=src_data.astype(np.float32),
                destination=destination,
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=dst_transform,
                dst_crs=target_crs,
                src_nodata=src_nodata,
                dst_nodata=src_nodata,
                resampling=Resampling[resampling_name],
            )

            valid_dst = destination != src_nodata

    rgba = apply_colormap(destination, valid_dst, data_min, data_max, stops)

    output_png_path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgba).save(output_png_path)
    print(f"Wrote PNG preview: {output_png_path} (using {sidecar_path.name})")


def build_composite_with_sidecar_and_preview(
    directory: Path,
    bbox: BoundingBox,
    start: YearMonth,
    end: YearMonth,
    tif_output_path: Path,
    sidecar_output_path: Path,
    png_output_path: Path,
    target_crs: CRS = WEBMERCATOR_CRS,
    resampling: str = "bilinear",
    zoom_resampling: str = "linear",
) -> None:
    """
    End-to-end pipeline: discover files, validate scaling, crop, average,
    write GeoTIFF composite, write colormap sidecar (including the real
    image_bounds and the client-side zoom_resampling hint), then render
    PNG preview.

    The crop window is computed once from the first file's own SDU transform
    and reused for every subsequent file, on the assumption that all monthly
    files in a run share the same grid (true for a single CM SAF product).
    """
    
    files = discover_monthly_files(directory, start, end)
    print(f"Using {len(files)} monthly files from {start.year}-{start.month:02d} to {end.year}-{end.month:02d}")

    scale_offset = validate_consistent_scaling(files)

    first_sdu_uri = resolve_sdu_uri(files[0])
    window = compute_crop_window(first_sdu_uri, bbox)
    edges = window_edges(first_sdu_uri, window)

    if (edges.west, edges.south, edges.east, edges.north) != (bbox.west, bbox.south, bbox.east, bbox.north):
        print(
            f"Note: requested bbox ({bbox.west}, {bbox.south}, {bbox.east}, {bbox.north}) was clamped/snapped "
            f"to the grid's real pixel edges ({edges.west}, {edges.south}, {edges.east}, {edges.north}). "
            "Downstream consumers (e.g. a web map) must use these real edges, not the requested bbox."
        )

    annual_grid, valid = compute_average_annual_sdu(files, window, scale_offset)
    if not valid.any():
        raise ValueError("Average annual crop contains no finite, non-nodata SDU values.")

    data_min, data_max = write_annual_grid_geotiff(
        annual_grid, valid, edges, files, start, end, scale_offset, tif_output_path
    )

    write_colormap_sidecar(
        tif_output_path, data_min, data_max, NODATA_VALUE, edges,
        target_crs, resampling, zoom_resampling, sidecar_output_path,
    )

    render_png_from_sidecar(tif_output_path, sidecar_output_path, png_output_path)


def main() -> None:
    directory = Path("./ORD69093").resolve()
    bbox = BoundingBox(west=-40.0, south=30.0, east=45.0, north=72.0)
    start = YearMonth(2021, 1)
    end = YearMonth(2026, 1)

    tif_output_path = Path("./output/sdu_europe_average_annual.tif")
    sidecar_output_path = Path("./output/sdu_europe_average_annual.colormap.json")
    png_output_path = Path("./output/sdu_europe_average_annual.png")

    build_composite_with_sidecar_and_preview(
        directory, bbox, start, end,
        tif_output_path, sidecar_output_path, png_output_path,
        target_crs=WEBMERCATOR_CRS, resampling="bilinear", zoom_resampling="linear",
    )


if __name__ == "__main__":
    main()
