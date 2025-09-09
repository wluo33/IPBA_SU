import arcpy
import os
import time
import math

# === Configuration ===
# Set environment
# change your path accordingly
# see if __name__ == "__main__" toward the bottom for additional input changes 
arcpy.env.overwriteOutput = True
arcpy.env.workspace = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt.gdb"
arcpy.env.scratchWorkspace = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt.gdb"

def clear_memory_workspace():
    original_ws = arcpy.env.workspace  # Backup current workspace
    arcpy.env.workspace = "in_memory"

    mem_fcs = arcpy.ListFeatureClasses() or []
    mem_rasters = arcpy.ListRasters() or []

    for name in mem_fcs + mem_rasters:
        try:
            arcpy.Delete_management(name)
        except:
            pass

    arcpy.env.workspace = original_ws  # Restore original workspace


def calculate_aspect_from_dem(dem_raster, aspect_raster):
    from arcpy.sa import Aspect
    aspect = Aspect(dem_raster)
    aspect.save(aspect_raster)
    return aspect_raster


def zonal_stats_aspect_dict(polygon_fc, aspect_raster, workspace, existing_stats=None, oids_to_update=None):
    out_table = os.path.join(workspace, "zonal_stats_tbl")
    if arcpy.Exists(out_table):
        arcpy.Delete_management(out_table)

    oid_field = "zone_id"
    if "zone_id" not in [f.name for f in arcpy.ListFields(polygon_fc)]:
        arcpy.management.AddField(polygon_fc, "zone_id", "LONG")
        with arcpy.da.UpdateCursor(polygon_fc, ["OBJECTID", "zone_id"]) as cursor:
            for oid, _ in cursor:
                cursor.updateRow([oid, oid])


    arcpy.sa.ZonalStatisticsAsTable(
        in_zone_data=polygon_fc,
        #zone_field=oid_field,
        zone_field="zone_id",
        in_value_raster=aspect_raster,
        out_table=out_table,
        ignore_nodata="DATA",
        statistics_type="ALL",
        circular_calculation="CIRCULAR",
        circular_wrap_value=360
    )

    stats = existing_stats.copy() if existing_stats else {}

    if "C_STD" not in [f.name for f in arcpy.ListFields(polygon_fc)]:
        arcpy.management.AddField(polygon_fc, "C_STD", "DOUBLE")
    if "C_MEAN" not in [f.name for f in arcpy.ListFields(polygon_fc)]:
        arcpy.management.AddField(polygon_fc, "C_MEAN", "DOUBLE")

    with arcpy.da.SearchCursor(out_table, [oid_field, "C_MEAN", "C_STD"]) as cursor:
        for oid, mean, std in cursor:
            if oids_to_update is None or oid in oids_to_update:
                stats[oid] = {"C_MEAN": mean, "C_STD": std}

    with arcpy.da.UpdateCursor(polygon_fc, [oid_field, "C_MEAN", "C_STD"]) as cursor:
        for row in cursor:
            oid = row[0]
            if oids_to_update is None or oid in oids_to_update:
                mean_val = stats.get(oid, {}).get("mean", None)
                std_val = stats.get(oid, {}).get("std", None)
                # Fallback to 0.0 if value is None
                row[1] = mean_val if mean_val is not None else 0.0
                row[2] = std_val if std_val is not None else 0.0
                cursor.updateRow(row)


    # Force ArcGIS to flush updated field values to schema (memory layers are picky)
    try:
        arcpy.management.CalculateField(polygon_fc, "C_MEAN", "!C_MEAN!", "PYTHON3")
        arcpy.management.CalculateField(polygon_fc, "C_STD", "!C_STD!", "PYTHON3")
    except Exception as e:
        print("⚠️ Warning while forcing field write:", e)



    print(f"✅ Assigned stats to {len(stats)} polygons")

    return stats



def get_neighbors(polygon_fc, workspace):
    neighbor_tbl = os.path.join(workspace, "neighbors")
    if arcpy.Exists(neighbor_tbl):
        arcpy.Delete_management(neighbor_tbl)

    arcpy.analysis.PolygonNeighbors(polygon_fc, neighbor_tbl, "", "NO_AREA_OVERLAP", "BOTH_SIDES")
    fields = [f.name for f in arcpy.ListFields(neighbor_tbl)]
    print("Neighbor table fields:", fields)
    try:
        src_field = next(f for f in fields if f.lower().startswith('src_'))
        nbr_field = next(f for f in fields if f.lower().startswith('nbr_'))
    except StopIteration:
        raise RuntimeError("Could not find suitable source and neighbor ID fields in PolygonNeighbors output.")

    neighbors = {}
    with arcpy.da.SearchCursor(neighbor_tbl, [src_field, nbr_field]) as cursor:
        for src, nbr in cursor:
            neighbors.setdefault(src, set()).add(nbr)
    return neighbors

def crosses_barrier(poly1, poly2, barrier_geoms, min_overlap_ratio=0.5):
    shared = poly1.boundary().intersect(poly2.boundary(), 2)
    if not shared or shared.length == 0:
        return False

    overlap_len = 0
    for barrier in barrier_geoms:
        inter = shared.intersect(barrier, 2)
        if inter and inter.length > 0:
            overlap_len += inter.length
    return (overlap_len / shared.length) > min_overlap_ratio



def merge_units(polygon_fc, aspect_raster, merge_below_area, max_aspect_diff, aspect_dict, min_area_pixel, ridge_geoms=None, valley_geoms=None, workspace="memory"):
    print("Building neighbor list...")
    neighbor_dict = get_neighbors(polygon_fc, workspace)

    print("Reading polygons...")
    geom_dict = {}
    area_dict = {}
    with arcpy.da.UpdateCursor(polygon_fc, ["OID@", "SHAPE@"]) as cursor:
        for oid, shape in cursor:
            geom_dict[oid] = shape
            area_dict[oid] = shape.area

    print("Sorting polygons by area ascending...")
    sorted_oids = sorted(area_dict.items(), key=lambda x: x[1])

    print("Identifying polygons for merging with global similarity sort...")
    scored_pairs = []
    used_ids = set()
    for oid, _ in sorted_oids:
        if oid not in geom_dict or oid in used_ids:
            continue
        if area_dict[oid] > merge_below_area:
            continue

        aspect = aspect_dict.get(oid, {}).get("C_MEAN")
        neighbors = neighbor_dict.get(oid, [])
        poly1 = geom_dict[oid]

        for nbr in neighbors:
            if nbr not in geom_dict or nbr in used_ids:
                continue

            poly2 = geom_dict[nbr]
            nbr_aspect = aspect_dict.get(nbr, {}).get("C_MEAN")
            if aspect is None or nbr_aspect is None:
                continue

            diff = min(abs(aspect - nbr_aspect), 360 - abs(aspect - nbr_aspect))
            if diff > max_aspect_diff:
                continue

            if ridge_geoms and crosses_barrier(poly1, poly2, ridge_geoms):
                continue
            if valley_geoms and crosses_barrier(poly1, poly2, valley_geoms):
                continue

            scored_pairs.append((diff, oid, nbr))




    scored_pairs.sort()
    merge_pairs = []
    for diff, oid, nbr in scored_pairs:
        if oid not in used_ids and nbr not in used_ids:
            merge_pairs.append((oid, nbr))
            used_ids.update([oid, nbr])

    print(f"Merging {len(merge_pairs)} pairs...")
    if not merge_pairs:
        return polygon_fc, set()  # ✅ Important fix: never return False

    # Assign merge IDs
    oid_to_merge_id = {}
    merge_id = 1
    for oid, nbr in merge_pairs:
        oid_to_merge_id[oid] = merge_id
        oid_to_merge_id[nbr] = merge_id
        merge_id += 1
    merged_oids = set(oid_to_merge_id.keys())

    merge_field = "merge_id"
    if merge_field in [f.name for f in arcpy.ListFields(polygon_fc)]:
        arcpy.management.DeleteField(polygon_fc, merge_field)
    arcpy.AddField_management(polygon_fc, merge_field, "LONG")

    with arcpy.da.UpdateCursor(polygon_fc, ["OID@", merge_field]) as cursor:
        for oid, _ in cursor:
            mid = oid_to_merge_id.get(oid)
            if mid:
                cursor.updateRow([oid, mid])

    # --- Merging ---
    out_fc = os.path.join(workspace, "merged_units")
    if arcpy.Exists(out_fc):
        arcpy.Delete_management(out_fc)

    print("Dissolving by merge_id...")
    merge_layer = os.path.join(workspace, "merge_layer")
    arcpy.MakeFeatureLayer_management(polygon_fc, merge_layer, f"{merge_field} IS NOT NULL")
    arcpy.Dissolve_management(merge_layer, out_fc, merge_field, multi_part="SINGLE_PART")

    print("Appending unmerged polygons...")
    merged_oids_flat = list(merged_oids)
    if merged_oids_flat:
        where_clause = f"NOT OBJECTID IN ({','.join(map(str, merged_oids_flat))})"
    else:
        where_clause = "1=1"  # all are unmerged

    unmerged_layer = os.path.join(workspace, "unmerged")
    arcpy.MakeFeatureLayer_management(polygon_fc, unmerged_layer, where_clause)
    arcpy.management.Append(unmerged_layer, out_fc, "NO_TEST")

    # --- Force singlepart & cleanup ---
    single_fc = os.path.join(workspace, "single_fc")
    arcpy.MultipartToSinglepart_management(out_fc, single_fc)
    arcpy.management.RepairGeometry(single_fc)

    eliminate_small_polygons(single_fc, aspect_raster, min_area_pixel=min_area_pixel)

    if merge_field in [f.name for f in arcpy.ListFields(single_fc)]:    
        arcpy.DeleteField_management(single_fc, merge_field)

    polygon_fc = arcpy.CopyFeatures_management(single_fc, os.path.join(workspace, "working_merged")).getOutput(0)

    if "zone_id" not in [f.name for f in arcpy.ListFields(polygon_fc)]:
        arcpy.management.AddField(polygon_fc, "zone_id", "LONG")

    with arcpy.da.UpdateCursor(polygon_fc, ["OBJECTID", "zone_id"]) as cursor:
        for oid, _ in cursor:
            cursor.updateRow([oid, oid])

    print(f"DEBUG: Polygon count after merge: {int(arcpy.management.GetCount(polygon_fc)[0])}")
    return polygon_fc, merged_oids





def compute_mean_circular_std_from_stats(aspect_stats):
    """Compute the mean circular std from the updated stats dictionary."""
    stds = [stats['C_STD'] for stats in aspect_stats.values() if stats['C_STD'] is not None]
    if stds:
        return sum(stds) / len(stds)
    else:
        return 0


def eliminate_small_polygons(fc, dem_raster, min_area_pixel, area_field="AREA"):
    desc = arcpy.Describe(dem_raster)
    real_min_area = min_area_pixel * desc.meanCellHeight * desc.meanCellWidth

    print(f"Minimum area threshold: {real_min_area:.2f} map units")

    working_fc = os.path.join("in_memory", "working_fc")
    arcpy.CopyFeatures_management(fc, working_fc)

    # Ensure area field exists
    existing_fields = [f.name for f in arcpy.ListFields(working_fc)]
    if area_field not in existing_fields:
        arcpy.AddField_management(working_fc, area_field, "DOUBLE")

    while True:
        # Recalculate area
        arcpy.management.CalculateField(working_fc, area_field, "!SHAPE.area!", "PYTHON3")

        # Select small polygons
        small_layer = "small_poly_lyr"
        where_clause = f"{area_field} <= {real_min_area}"
        arcpy.MakeFeatureLayer_management(working_fc, small_layer, where_clause)
        count = int(arcpy.GetCount_management(small_layer)[0])

        if count == 0:
            print("No polygons smaller than threshold. Finished elimination.")
            break

        print(f"Eliminating {count} polygons smaller than threshold...")
        eliminated_fc = os.path.join("in_memory", "eliminated")
        arcpy.management.Eliminate(small_layer, eliminated_fc, "LENGTH")

        # Convert to singlepart
        singlepart_fc = os.path.join("in_memory", "singlepart")
        arcpy.MultipartToSinglepart_management(eliminated_fc, singlepart_fc)

        # Replace working_fc with updated version
        arcpy.Delete_management(working_fc)
        arcpy.CopyFeatures_management(singlepart_fc, working_fc)

        # Clean temp
        for temp in [small_layer, eliminated_fc, singlepart_fc]:
            if arcpy.Exists(temp):
                arcpy.Delete_management(temp)

    # Final overwrite to fc
    arcpy.Delete_management(fc)
    arcpy.CopyFeatures_management(working_fc, fc)
    arcpy.Delete_management(working_fc)
    print(f"Final cleaned polygons saved to: {fc}")


def fix_enclosed_islands(polygon_fc, workspace):
    print("🔁 Fixing enclosed island polygons...")
    if "enclose_id" in [f.name for f in arcpy.ListFields(polygon_fc)]:
        arcpy.DeleteField_management(polygon_fc, "enclose_id")
    arcpy.AddField_management(polygon_fc, "enclose_id", "LONG")

    polygons = []
    with arcpy.da.SearchCursor(polygon_fc, ["OID@", "SHAPE@"]) as cursor:
        for row in cursor:
            polygons.append((row[0], row[1]))

    with arcpy.da.UpdateCursor(polygon_fc, ["OID@", "SHAPE@", "enclose_id"]) as cursor:
        for i, (oid_a, shape_a, _) in enumerate(cursor):
            for oid_b, shape_b in polygons:
                if oid_a == oid_b:
                    continue
                if shape_a.within(shape_b):
                    cursor.updateRow([oid_a, shape_a, oid_b])
                    break



    enclosed_count = 0
    with arcpy.da.SearchCursor(polygon_fc, ["enclose_id"]) as cursor:
        enclosed_count = sum(1 for row in cursor if row[0] is not None)

    if enclosed_count == 0:
        print("⚠️ No enclosed polygons found — skipping dissolve.")
        return polygon_fc

    arcpy.Dissolve_management(polygon_fc, enclosed_fc, "enclose_id", multi_part="SINGLE_PART")


    remaining_layer = arcpy.management.MakeFeatureLayer(polygon_fc, "remaining_layer", "enclose_id IS NULL")
    arcpy.Append_management(remaining_layer, enclosed_fc, "NO_TEST")

    cleaned_fc = os.path.join(workspace, "merged_after_enclosure_fix")
    return arcpy.CopyFeatures_management(enclosed_fc, cleaned_fc).getOutput(0)


def update_polygon_fields(polygon_fc, aspect_stats):
    """Update C_MEAN and C_STD fields using aspect_stats dictionary."""
    fields = ["OBJECTID", "C_MEAN", "C_STD"]

    # Ensure fields exist
    existing_fields = [f.name for f in arcpy.ListFields(polygon_fc)]
    if "C_MEAN" not in existing_fields:
        arcpy.AddField_management(polygon_fc, "C_MEAN", "DOUBLE")
    if "C_STD" not in existing_fields:
        arcpy.AddField_management(polygon_fc, "C_STD", "DOUBLE")

    with arcpy.da.UpdateCursor(polygon_fc, fields) as cursor:
        for row in cursor:
            oid = row[0]
            if oid in aspect_stats:
                row[1] = aspect_stats[oid]['C_MEAN']
                row[2] = aspect_stats[oid]['C_STD']
                cursor.updateRow(row)



def iterative_merge(polygon_fc, aspect_raster, merge_below_area, max_aspect_diff, max_iters=10, std_threshold=45, min_area_pixel=15, ridge_geoms=None, valley_geoms=None, output_fc=None):
    print("Starting iterative merging...")
    workspace = "memory"
    polygon_fc = arcpy.management.CopyFeatures(polygon_fc, os.path.join(workspace, "working_fc")).getOutput(0)

    # 🔹 Compute initial aspect stats
    aspect_stats = zonal_stats_aspect_dict(polygon_fc, aspect_raster, workspace)
    update_polygon_fields(polygon_fc, aspect_stats)
    mean_std = compute_mean_circular_std_from_stats(aspect_stats)
    print(f"Initial mean circular std = {mean_std:.2f}")


    if mean_std > std_threshold:
        print(f"Initial circular std ({mean_std:.2f}) exceeds threshold ({std_threshold}). Stopping.")
        arcpy.CopyFeatures_management(polygon_fc, output_fc)
        return

    prev_mean_std = mean_std

    for i in range(max_iters):
        print(f"\n--- Iteration {i+1} ---")

        # 🔁 Merge units
        polygon_fc, merged_oids = merge_units(
            polygon_fc,
            aspect_raster,
            merge_below_area,
            max_aspect_diff,
            aspect_stats,
            min_area_pixel, 
            ridge_geoms=ridge_geoms,
            valley_geoms=valley_geoms,
            workspace=workspace
        )

        print("\nDEBUG: Polygon count after merge:", arcpy.GetCount_management(polygon_fc)[0])

        # 🔁 Optional fix of enclosed islands
        polygon_fc = fix_enclosed_islands(polygon_fc, workspace)

        if not merged_oids:
            print("No polygons merged this iteration. Stopping.")
            break

        # 🔁 Recompute aspect stats after merging
        aspect_stats = zonal_stats_aspect_dict(polygon_fc, aspect_raster, workspace)
        mean_std = compute_mean_circular_std_from_stats(aspect_stats)
        print(f"DEBUG: Circular stds this round: {[round(s['C_STD'],2) for s in aspect_stats.values() if s['C_STD'] is not None][:10]}...")


        curr_poly_count = int(arcpy.management.GetCount(polygon_fc)[0])

        print(f"Mean circular std = {mean_std:.2f}")
        print(f"Polygon count: {curr_poly_count}")
        print(f"Merged polygons this iteration: {len(merged_oids)}")

        if mean_std > std_threshold:
            print("Circular std exceeded threshold. Stopping.")
            break

        if abs(prev_mean_std - mean_std) < 1e-2:
            print("Mean circular std did not change significantly. Stopping.")
            break

        prev_mean_std = mean_std


        # 🛠 Optional debug export
        debug_fc = f"D:\\ibpa_slpu\\gpt_new\\paper2025\\new_clean_tpi\\debug_iter_{i+1}.shp"
        update_polygon_fields(polygon_fc, aspect_stats)  # 🧠 Add this line before saving
        arcpy.CopyFeatures_management(polygon_fc, debug_fc)



        # 🧹 Clear memory
        clear_memory_workspace()

        # 🧾 Area summary
        print("💡 Unique area values after merge:")
        with arcpy.da.SearchCursor(polygon_fc, ["SHAPE@AREA"]) as cursor:
            areas = sorted([row[0] for row in cursor])
            print(f" - Min area: {areas[0]:.2f}, Max area: {areas[-1]:.2f}, Count: {len(areas)}")




    # 🧼 Final cleanup
    print("Final cleanup: eliminating small polygons...")
    eliminate_small_polygons(polygon_fc, aspect_raster, min_area_pixel=25)

    # 🔁 Final aspect stats
    print("🔄 Recomputing final C_MEAN and C_STD before saving...")
    aspect_stats = zonal_stats_aspect_dict(polygon_fc, aspect_raster, workspace)
    update_polygon_fields(polygon_fc, aspect_stats)
    final_mean_std = compute_mean_circular_std_from_stats(aspect_stats)
    print(f"Final mean circular std = {final_mean_std:.2f}")



    # ✨ Write C_MEAN and C_STD to attribute table
    with arcpy.da.UpdateCursor(polygon_fc, ["OBJECTID", "C_MEAN", "C_STD"]) as cursor:
        for row in cursor:
            oid = row[0]
            if oid in aspect_stats:
                row[1] = aspect_stats[oid]["C_MEAN"]
                row[2] = aspect_stats[oid]["C_STD"]
                cursor.updateRow(row)


    # 💾 Save result
    if output_fc:
        print(f"Saving final output to: {output_fc}")
        arcpy.CopyFeatures_management(polygon_fc, output_fc)

        # ✅ Optional: Check if output has non-zero values
        with arcpy.da.SearchCursor(output_fc, ["C_STD"]) as cursor:
            std_vals = [row[0] for row in cursor if row[0] is not None]
            if all(s == 0 for s in std_vals):
                print("⚠️ All C_STD values in saved output are zero!")
            else:
                print("✅ Final output contains valid C_STD values.")

    print("✅ Iterative merging complete.")





if __name__ == "__main__":
    start = time.time()

    #additional input/output variables
    #modify to reflect your own cases
    slope_units = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt_all1.shp"                #ouput shapefile from IPBA
    dem_raster = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt.tif"                      #input DEM
    output_fc = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt_all1_75k_30_40_30.shp"     #output merged slope unit shapefile


    merge_below_area = 75000            #in square meter (m^2), polygons with area below this will be considered for merge
    max_aspect_diff = 30                #aspect difference below this will considered for merge
    max_iters = 10                      #maximum number of iterations
    std_threshold=40                    #overall standard deviation of aspect threshold to stop iteration
    min_area_pixel = 30                 #polygon smaller than this number of pixels will be merged with neighbor regardless of aspect


    
    rigval_buf = 30 #default to pixel size,
                    #here the pixel size is 30 m, so a buffer of 30 m is used to test if the border between 2 neighboring polygons is a ridge/valley

##  since we are using the accentuated DEM, the following section is unnecessary; if you want to check if the merging is crossing ridge/valley line barrier,
##  you may uncomment the following section.
    
##    ridge_raster = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt_ridge_rg.tif"
##    valley_raster = r"D:\ibpa_slpu\paper\code_for_submit\belm\bt_valley_rg.tif"
##    ridge_lines = "memory\\ridges"
##    valley_lines = "memory\\valleys"
##    arcpy.conversion.RasterToPolyline(ridge_raster, ridge_lines, "ZERO", 0, "NO_SIMPLIFY")
##    arcpy.conversion.RasterToPolyline(valley_raster, valley_lines, "ZERO", 0, "NO_SIMPLIFY")
##
##    ridge_polygons = "memory\\ridges_poly"
##    valley_polygons = "memory\\valleys_poly"
##
##    arcpy.analysis.Buffer(ridge_lines, ridge_polygons, rigval_buf)
##    arcpy.analysis.Buffer(valley_lines, valley_polygons, rigval_buf)
##    
##    ridge_geoms = [row[0] for row in arcpy.da.SearchCursor(ridge_polygons, ["SHAPE@"])]
##    valley_geoms = [row[0] for row in arcpy.da.SearchCursor(valley_polygons, ["SHAPE@"])]



    aspect_raster = r"memory\aspect"    #hold aspect raster, but store in memory


    # Copy input to avoid modifying original
    input_copy = r"memory\input_copy"
    arcpy.management.CopyFeatures(slope_units, input_copy)
    input_fc = input_copy
    
   
    print("Computing aspect from DEM...")
    calculate_aspect_from_dem(dem_raster, aspect_raster)
    

##call iterative merge; by default, there is no cross ridge/valley check (ridge_geoms=None, valley_geoms=None),
##because using accentuate DEM making cross check unnecessary.
##if you uncomment out the section mentioned above, you should use the following call instead
##    iterative_merge(
##        polygon_fc=input_fc,
##        aspect_raster=aspect_raster,
##        merge_below_area=merge_below_area,
##        max_aspect_diff=max_aspect_diff,
##        max_iters=max_iters,
##        std_threshold=std_threshold,
##        min_area_pixel=min_area_pixel,
##        ridge_geoms = ridge_geoms,
##        valley_geoms = valley_geoms,
##        output_fc=output_fc       
##    )    
##def iterative_merge(polygon_fc, aspect_raster, merge_below_area, max_aspect_diff, max_iters=10, std_threshold=45, min_area_pixel=15, ridge_geoms=None, valley_geoms=None, output_fc=None):
    iterative_merge(
        polygon_fc=input_fc,
        aspect_raster=aspect_raster,
        merge_below_area=merge_below_area,
        max_aspect_diff=max_aspect_diff,
        max_iters=max_iters,
        std_threshold=std_threshold,
        min_area_pixel=min_area_pixel,
        output_fc=output_fc
    )

    print(f"\n\U0001F552 Total processing time: {round(time.time() - start, 2)} seconds")


