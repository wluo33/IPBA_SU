import arcpy
from arcpy.sa import *
import os

import networkx as nx

#need to install NetworkX first (e.g., conda install -c conda-forge networkx)


# === Configuration ===
# Set environment
# change your path accordingly
output_path = r"D:\ibpa_slpu\paper\code_for_submit\belm"
workspace = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt.gdb"
input_dem0 = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt.tif"


L2 = 15                     #smoothing window size (ideally odd number)
L = 15                      #search radius for geomorphon
                            #now set both L and L2 the same by default, user may experiment, see discussion in the paper

MIN_PIXEL_COUNT = 25        # Minimum number of pixels for valid polygon parts



# Input DEM
input_dem = Raster(input_dem0)


scratch_workspace = workspace
arcpy.env.workspace = workspace
arcpy.env.scratchWorkspace = scratch_workspace
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")



base_name0 = os.path.splitext(os.path.basename(input_dem0))[0]
print("base_name0=",base_name0)

# Area/size thresholds
cell_x = float(arcpy.GetRasterProperties_management(input_dem, "CELLSIZEX").getOutput(0))
cell_y = float(arcpy.GetRasterProperties_management(input_dem, "CELLSIZEY").getOutput(0))
  
#HOLE_THRESHOLD = MIN_PIXEL_COUNT * cell_x * cell_y        # In map units^2 (e.g., square meters)
HOLE_THRESHOLD = MIN_PIXEL_COUNT * cell_x * cell_y        # In map units^2 (e.g., square meters)
MIN_PRUNE_LENGTH = MIN_PIXEL_COUNT * cell_x

# Output paths
smooth_dem = "smooth_dem"
geomorphon_output = os.path.join(output_path, "geomorphon_output.tif")
valley_output = os.path.join(output_path, f"{base_name0}_valley_output.tif")
ridge_output = os.path.join(output_path, f"{base_name0}_ridge_output.tif")


# Neighborhood settings (e.g., 15x15 window)
neighborhood = NbrRectangle(L2, L2, "CELL")

# Compute mean elevation of neighborhood
mean_elev = FocalStatistics(input_dem, neighborhood, "MEAN", "DATA")

# TPI = Cell Elevation - Mean Neighborhood Elevation
tpi = input_dem + 10 * (input_dem - mean_elev)
#tpi = input_dem + (input_dem - mean_elev)
tpi.save("tpi")

# === Preprocessing: Smooth TPI ===
neighborhood2 = NbrRectangle(7, 7, "CELL")
smooth_dem = FocalStatistics(tpi, neighborhood2, "MEAN", "DATA")
smooth_dem.save("smooth_dem")

# === Geomorphon Analysis ===
print("Running GeomorphonLandforms...")
geomorphon = GeomorphonLandforms(smooth_dem, "", "1", "CELLS", L, "", "")
geomorphon.save(geomorphon_output)

valley_output0 = "valley_output0"
# Extract valleys (pit = 9, valley = 10)
arcpy.gp.Con_sa(geomorphon_output, "1", valley_output0, "", "VALUE = 9 Or VALUE = 10")
valley_output = FocalStatistics(valley_output0, NbrRectangle(3, 3, "CELL"), "MAJORITY", "DATA")

ridge_output0 = "ridge_output0"
# Extract ridges (ridge = 2, peak = 3)
arcpy.gp.Con_sa(geomorphon_output, "1", ridge_output0, "", "VALUE = 2 Or VALUE = 3")
ridge_output = FocalStatistics(ridge_output0, NbrRectangle(3, 3, "CELL"), "MAJORITY", "DATA")


def remove_isolated_lines_precise(input_fc, output_fc, pixel_size=None, pixel_threshold=15):
    """
    Removes isolated line features from input_fc.
    Keeps isolated lines if they are longer than pixel_threshold * pixel_size.
    Saves connected lines and long isolated lines to output_fc.
    """

    # Step 1: Determine pixel size if not provided
    if pixel_size is None:
        desc = arcpy.Describe(input_fc)
        if hasattr(desc, "extent") and hasattr(desc, "spatialReference"):
            cell_x = float(arcpy.GetRasterProperties_management("b.tif", "CELLSIZEX").getOutput(0))
            cell_y = float(arcpy.GetRasterProperties_management("b.tif", "CELLSIZEY").getOutput(0))
            pixel_size = max(cell_x, cell_y)
        else:
            raise ValueError("Raster reference required to calculate pixel size.")

    length_threshold = pixel_threshold * pixel_size

    # Step 2: Add a temporary unique ID if needed
    temp_id_field = "TEMP_ID"
    if temp_id_field not in [f.name for f in arcpy.ListFields(input_fc)]:
        arcpy.management.AddField(input_fc, temp_id_field, "LONG")
        with arcpy.da.UpdateCursor(input_fc, ["OID@", temp_id_field]) as cursor:
            for oid, _ in cursor:
                cursor.updateRow([oid, oid])

    # Step 3: Use Feature Vertices To Points to get endpoints
    endpoints = os.path.join("in_memory", "endpoints")
    arcpy.management.FeatureVerticesToPoints(input_fc, endpoints, "BOTH_ENDS")

    # Step 4: Count how many endpoints each line shares
    joined = os.path.join("in_memory", "joined")
    arcpy.analysis.SpatialJoin(endpoints, input_fc, joined, "JOIN_ONE_TO_MANY", match_option="INTERSECT")

    # Step 5: Count number of joins per line
    endpoint_counts = {}
    with arcpy.da.SearchCursor(joined, [temp_id_field]) as cursor:
        for (temp_id,) in cursor:
            endpoint_counts[temp_id] = endpoint_counts.get(temp_id, 0) + 1

    # Step 6: Select connected lines
    connected_ids = {tid for tid, count in endpoint_counts.items() if count > 2}

    # Step 7: Create final selection (connected + isolated longer than threshold)
    arcpy.MakeFeatureLayer_management(input_fc, "lines_lyr")

    long_isolated_ids = []
    with arcpy.da.SearchCursor(input_fc, [temp_id_field, "SHAPE@LENGTH"]) as cursor:
        for tid, length in cursor:
            if tid not in connected_ids and length > length_threshold:
                long_isolated_ids.append(tid)

    all_keep_ids = connected_ids.union(long_isolated_ids)

    where_clause = f"{temp_id_field} IN ({','.join(map(str, all_keep_ids))})"
    arcpy.management.SelectLayerByAttribute("lines_lyr", "NEW_SELECTION", where_clause)

    arcpy.management.CopyFeatures("lines_lyr", output_fc)
    print(f"Saved output with {len(all_keep_ids)} features to: {output_fc}")





def prune_and_merge_trunks(line_fc, output_fc, length_threshold=1000):
    """
    Prune short branches and merge remaining trunk segments into continuous lines.
    
    Parameters:
    - line_fc (str): Input line shapefile or feature class.
    - output_fc (str): Output merged shapefile after pruning short branches.
    - length_threshold (float): Branches shorter than this will be pruned.
    """
    G = nx.Graph()
    line_dict = {}

    # Build graph and store line geometries
    with arcpy.da.SearchCursor(line_fc, ['OID@', 'SHAPE@']) as cursor:
        for oid, geom in cursor:
            start = (round(geom.firstPoint.X, 3), round(geom.firstPoint.Y, 3))
            end = (round(geom.lastPoint.X, 3), round(geom.lastPoint.Y, 3))
            G.add_edge(start, end, oid=oid, length=geom.length)
            line_dict[oid] = geom

    # Identify and remove short branches
    endpoints = [n for n in G.nodes if G.degree[n] == 1]
    oids_to_remove = set()
    for ep in endpoints:
        path = [ep]
        total_length = 0
        current = ep
        while True:
            neighbors = list(G.neighbors(current))
            if len(neighbors) != 1:
                break
            next_node = neighbors[0]
            edge_data = G.get_edge_data(current, next_node)
            total_length += edge_data['length']
            path.append(next_node)
            current = next_node
            if G.degree[next_node] != 2:
                break
        if total_length < length_threshold:
            for i in range(len(path) - 1):
                oid = G[path[i]][path[i + 1]]['oid']
                oids_to_remove.add(oid)

    # Remove short branches from graph
    G.remove_edges_from([
        (u, v) for u, v, d in G.edges(data=True) if d['oid'] in oids_to_remove
    ])

    # Traverse and merge remaining trunk segments
    merged_geoms = []
    visited_edges = set()

    for u in G.nodes:
        if G.degree[u] != 2:
            for v in G.neighbors(u):
                edge = (u, v)
                if edge in visited_edges or (v, u) in visited_edges:
                    continue
                path = [u, v]
                current = v
                prev = u
                while G.degree[current] == 2:
                    next_nodes = list(G.neighbors(current))
                    next_nodes.remove(prev)
                    prev = current
                    current = next_nodes[0]
                    path.append(current)
                    if G.degree[current] != 2:
                        break
                # Rebuild merged geometry
                merged = line_dict[G[path[0]][path[1]]['oid']]
                for i in range(1, len(path) - 1):
                    oid = G[path[i]][path[i + 1]]['oid']
                    merged = merged.union(line_dict[oid])
                    visited_edges.add((path[i], path[i + 1]))
                visited_edges.add((path[0], path[1]))
                merged_geoms.append(merged)


    # Save to output_fc
    # Create empty output feature class
    if not arcpy.Exists("in_memory\\empty_line"):
        spatial_ref = arcpy.Describe(line_fc).spatialReference
        arcpy.CreateFeatureclass_management("in_memory", "empty_line", "POLYLINE", spatial_reference=spatial_ref)

    arcpy.CopyFeatures_management("in_memory\\empty_line", output_fc)

    with arcpy.da.InsertCursor(output_fc, ['SHAPE@']) as ic:
        for geom in merged_geoms:
            ic.insertRow([geom])





# === Function: Process Ridge/Valley Mask ===
def process_valley_or_ridge(mask_raster, out_name, hole_thresh=HOLE_THRESHOLD, min_pixels=MIN_PIXEL_COUNT):
    print(f"Processing: {out_name}")

 


    # Step 1: Raster to Polygon
    mask_polygon = f"{out_name}_poly"
    with arcpy.EnvManager(outputZFlag="Disabled", outputMFlag="Disabled"):
        arcpy.conversion.RasterToPolygon(
            in_raster=mask_raster,
            out_polygon_features=mask_polygon,
            simplify="NO_SIMPLIFY",
            raster_field="Value",
            create_multipart_features="SINGLE_OUTER_PART"
        )

    # Step 2: Eliminate small internal holes
    mask_poly_elim = f"{out_name}_poly_e"
    arcpy.management.EliminatePolygonPart(
        in_features=mask_polygon,
        out_feature_class=mask_poly_elim,
        condition="AREA",
        part_area=hole_thresh,
        part_area_percent=0,
        part_option="ANY"
    )

    # Step 3: Remove small polygons
    area_thresh = HOLE_THRESHOLD
    
    singleparts_fc = f"{out_name}_singleparts"
    arcpy.MultipartToSinglepart_management(mask_poly_elim, singleparts_fc)

    lyr = "layer"
    arcpy.management.MakeFeatureLayer(singleparts_fc, lyr)
    arcpy.management.SelectLayerByAttribute(lyr, "NEW_SELECTION", f"Shape_Area < {area_thresh}")
    arcpy.management.DeleteFeatures(lyr)

    # Step 4: Dissolve to simplify topology and save result
    dissolved = f"{out_name}_polygon_cleaned"
    arcpy.management.Dissolve(lyr, dissolved)

    
    
    print(f"✔ Cleaned polygon features saved to: {dissolved}")


    # Step 5: Polygon to Centerline
    centerline = f"{out_name}_centerline" 

    arcpy.topographic.PolygonToCenterline(dissolved, centerline)
    arcpy.CopyFeatures_management(centerline, os.path.join(output_path, f"{out_name}_centerline.shp"))


    centerline_pruned = f"{out_name}_centerline_pruned"
    centerline_raster = f"{out_name}_centerline_raster"
    min_valrig_length = MIN_PRUNE_LENGTH
    prune_and_merge_trunks(centerline, centerline_pruned, min_valrig_length)
    #def prune_and_merge_trunks(line_fc, output_fc, length_threshold=1000):

    #prune_short_branches(centerline, centerline_pruned, min_valrig_length)
    #def prune_short_branches(line_fc, output_fc, length_threshold=1000):

    with arcpy.EnvManager(snapRaster=input_dem, extent=input_dem):
        arcpy.conversion.FeatureToRaster(
            in_features=centerline_pruned,
            field="OBJECTID",
            out_raster=centerline_raster,
            cell_size=input_dem,
        )

     
    #center_raster2 = Raster(centerline_raster)
    Raster(centerline_raster).save(out_name)

    print(f"✔ Final output raster: {out_name}")



# === Process Valley ===
      
process_valley_or_ridge(valley_output, "valley_rg")

arcpy.management.CopyRaster(
    "valley_rg",
    os.path.join(output_path, f"{base_name0}_valley_rg.tif"),
    pixel_type = "32 bit signed",
    nodata_value="2147483647",
    format="TIFF"
)




# === Process Ridge ===

process_valley_or_ridge(ridge_output, "ridge_rg")

arcpy.management.CopyRaster(
    "ridge_rg",
    os.path.join(output_path, f"{base_name0}_ridge_rg.tif"),
    pixel_type = "32 bit signed",
    nodata_value="2147483647",
    format="TIFF"
)



print("✅ Done processing all valley and ridge rasters.")
