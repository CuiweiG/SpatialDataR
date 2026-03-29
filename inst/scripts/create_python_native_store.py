#!/usr/bin/env python3
"""Create a Python-native SpatialData Zarr store from MERFISH data.

This produces a store identical to what spatialdata-io would generate,
using Parquet for points and proper Zarr metadata.

Data: Allen MERFISH VISp (Moffitt et al. 2018 Science)
"""
import pandas as pd
import numpy as np
from pathlib import Path
import spatialdata as sd
from spatialdata.models import PointsModel, TableModel, ShapesModel
import geopandas as gpd
from shapely.geometry import Point
import anndata as ad

print("spatialdata version:", sd.__version__)

# Read raw MERFISH data
csv_path = Path("C:/Users/win10/merfish_real.csv")
print(f"Reading {csv_path} ...")
df = pd.read_csv(csv_path)
print(f"  {len(df)} rows, columns: {list(df.columns)}")

# Subset to bounding box (same as spatialdata-sandbox)
bb = {"x0": 1154, "x1": 3172, "y0": 4548, "y1": 6566}
sub = df[(df.x_um > bb["x0"]) & (df.x_um < bb["x1"]) &
         (df.y_um > bb["y0"]) & (df.y_um < bb["y1"])].copy()
print(f"  Subsetted: {len(sub)} spots")

# Subsample for manageable size (100K)
if len(sub) > 100000:
    sub = sub.sample(100000, random_state=42)
    print(f"  Subsampled to {len(sub)}")

# Create Points element (Parquet-native)
points_df = pd.DataFrame({
    "x": sub.x_um.values,
    "y": sub.y_um.values,
    "gene": pd.Categorical(sub.gene.values),
})
points = PointsModel.parse(
    points_df,
    coordinates={"x": "x", "y": "y"},
    feature_key="gene",
)
print(f"  Points: {len(points)} transcripts")

# Create Shapes (cell boundaries as circles)
layers = sub.layer.unique()
cells = []
cell_id = 0
for layer in layers:
    lsub = sub[sub.layer == layer]
    # Sample 20 representative positions per layer
    n = min(20, len(lsub))
    sample = lsub.sample(n, random_state=42)
    for _, row in sample.iterrows():
        cells.append({
            "geometry": Point(row.x_um, row.y_um),
            "radius": 10.0,
            "cell_id": cell_id,
            "cell_type": layer,
        })
        cell_id += 1

gdf = gpd.GeoDataFrame(cells)
shapes = ShapesModel.parse(gdf, transformations={})
print(f"  Shapes: {len(shapes)} cells, {len(layers)} layers")

# Create Table (cell metadata as AnnData)
obs = pd.DataFrame({
    "cell_type": [c["cell_type"] for c in cells],
    "cell_id": [c["cell_id"] for c in cells],
})
obs.index = obs.index.astype(str)
adata = ad.AnnData(obs=obs)
adata.obs["region"] = "cell_boundaries"
adata.obs["instance_key"] = range(len(cells))
table = TableModel.parse(
    adata,
    region="cell_boundaries",
    region_key="region",
    instance_key="instance_key",
)
print(f"  Table: {adata.n_obs} cells")

# Build SpatialData object
sdata = sd.SpatialData(
    points={"transcripts": points},
    shapes={"cell_boundaries": shapes},
    tables={"table": table},
)
print(f"\nSpatialData object created successfully")

# Write to Zarr
out_path = "C:/Users/win10/merfish_python_native.zarr"
print(f"\nWriting to {out_path} ...")
sdata.write(out_path, overwrite=True)

# Verify
import zarr
store = zarr.open(out_path, mode="r")
print(f"\nZarr store contents:")
def print_tree(group, prefix=""):
    for key in sorted(group.keys()):
        item = group[key]
        if hasattr(item, 'keys'):
            print(f"  {prefix}{key}/")
            if len(prefix) < 8:
                print_tree(item, prefix + "  ")
        else:
            print(f"  {prefix}{key}: {item.shape} {item.dtype}")

try:
    print_tree(store)
except:
    print("  (tree print skipped due to encoding)")
print(f"\nDone. Store at: {out_path}")
