import splusdata
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

conn = splusdata.Core()

bands = ["r", "g", "i", "J0378", "J0395", "J0410", "J0430", "g", "J0515", "r", "J0660", "i", "J0861"]

tables = [
    "idr4_dual.idr4_detection_image",
    # "idr4_dual.idr4_dual_%",
    "idr4_single.idr4_single_%",
    "idr4_psf.idr4_psf_%"
]

def get_tables():
    tabs = []
    for table in tables:
        if '%' in table:
            for band in bands:
                tabs.append(table.replace("%", band))
        else:
            tabs.append(table)
    return tabs

def format_ra_dec(ra, dec):
    # Create a SkyCoord object
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    # Format RA and DEC in the desired format
    # RA: hours, minutes, seconds; DEC: degrees, arcminutes, arcseconds
    # 'precision' controls the number of decimal places
    # 'pad=True' ensures leading zeros are included
    # 'alwayssign=True' for DEC ensures the + sign is included for positive values
    ra_formatted = coord.ra.to_string(unit=u.hour, sep='', precision=2, pad=True)
    dec_formatted = coord.dec.to_string(sep='', precision=1, alwayssign=True, pad=True)

    return f"SPLUS{ra_formatted}{dec_formatted}"

def match_coordinates_parallel(df, overlap_df=None):
    df["ID"] = df["ID"].apply(lambda x: x.decode('utf-8'))
    IDs = df["ID"].tolist()

    if overlap_df is not None:
        combined_df = pd.concat([overlap_df, df]).drop_duplicates(subset=['ID']).reset_index(drop=True)
    else:
        combined_df = df

    ra_col, dec_col = get_ra_dec_columns(combined_df)
    combined_df['unique_id'] = combined_df.apply(lambda row: format_ra_dec(row[ra_col], row[dec_col]), axis=1)
    
    # Ensure the columns are numeric and not containing NaN values
    combined_df[ra_col] = pd.to_numeric(combined_df[ra_col], errors='coerce')
    combined_df[dec_col] = pd.to_numeric(combined_df[dec_col], errors='coerce')
    combined_df.dropna(subset=[ra_col, dec_col], inplace=True)
    
    combined_df.reset_index(drop=True, inplace=True)

    # Convert RA and DEC to astropy Quantity
    ra_quantity = combined_df[ra_col].values * u.degree
    dec_quantity = combined_df[dec_col].values * u.degree
    coords = SkyCoord(ra=ra_quantity, dec=dec_quantity)

    def process_row(i):
        ## Check for CLASS_STAR of object
        class_star_col = get_column_from_pattern(combined_df, "class_star")
        if combined_df.loc[i, class_star_col] > 0.5:
            match_radius = 5 * u.arcsec
        else:
            match_radius = 3 * u.arcsec

    # Adjust based on your parameter
        idx, d2d, _ = coords[i].match_to_catalog_sky(coords, nthneighbor=2)
        if d2d < match_radius and idx != i:
            return i, int(idx)

    # Parallel processing
    with ThreadPoolExecutor(max_workers=128) as executor:
        futures = {executor.submit(process_row, i): i for i in range(len(combined_df))}
        for future in as_completed(futures):
            result = future.result()
            if result:
                i, idx = result
                combined_df.iloc[idx, combined_df.columns.get_loc('unique_id')] = combined_df.iloc[i, combined_df.columns.get_loc('unique_id')]
    
    # Select only IDs from list IDs in another way 

    return combined_df[combined_df["ID"].isin(IDs)]

def get_ra_dec_columns(table):
    ra = None
    dec = None
    for column in table.columns:
        if 'ra' in column.lower():
            ra = column
        if 'dec' in column.lower():
            dec = column
    
    if ra is None or dec is None:
        raise Exception("Could not find RA/Dec columns in table.")
    return ra, dec

def get_column_from_pattern(table, pattern):
    for column in table.columns:
        if pattern in column.lower():
            return column
        
    raise Exception("Could not find column with pattern: ", pattern)

def create_query(table, ra=None, width=None):
    if 'psf' in table or 'single' in table:
        band = table.split("_")[-1]
        query = f"""
            SELECT id, id_ra as ra, id_dec as dec, class_star_{band} as class_star FROM {table}
            WHERE 1 = CONTAINS(POINT('ICRS', id_ra, id_dec), BOX('ICRS', {ra}, 0, {width}, 300))
        """
        
    else:
        query = f"""
            SELECT id, ra, dec, class_star FROM {table}
            WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), BOX('ICRS', {ra}, 0, {width}, 300))
        """
        
    res = conn.query(query).to_pandas()
    return res

def count_query(table, ra=None, width=None):
    query = f"""
           SELECT COUNT(id) FROM {table}
           WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), BOX('ICRS', {ra}, 0, {width}, 300))
           """
           
    res = conn.query(query)
    return res["COUNT"][0]



@dataclass 
class TableClass:
    tables = []

    def get_table(self, table, ra, width):
        print("Getting table: ", table)
        table = create_query(table, ra, width)
        self.tables.append(table)
    
    def get_concat_table(self):
        #remove Nan values
        for table in self.tables:
            table.dropna(inplace=True)

        return pd.concat(self.tables)

    def parallel_match_data(self, previous_df):
        return match_coordinates_parallel(self.get_concat_table(), previous_df)

    def reset_tables(self):
        self.tables = []



def parallel_get_data(data_obj, ra, width):
    thread_pool = ThreadPoolExecutor(max_workers=20)
    for table in get_tables():
        thread_pool.submit(data_obj.get_table, table, ra, width)
    thread_pool.shutdown(wait=True)

def parallel_get_data_wt_wait(thread_pool, tmp_data_obj, ra, width):
    for table in get_tables():
        thread_pool.submit(tmp_data_obj.get_table, table, ra, width)

def main():
    tmp_thread_pool = None
    advance = False
    tmp_data_obj = None

    ra_range = [0, 1]  # Adjust as needed
    ra = ra_range[0]
    width = 0.1  # Size of each query box
    overlap = 0.1

    previous_df = pd.DataFrame()
    data_obj = TableClass()
    
    print("Starting...")
    print("Creating Unique ID for tables: ", get_tables())

    ## Check for arguments --ra and --width
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--initra", help="Initial RA in degrees", type=float)
    parser.add_argument("--width", help="Width in degrees", type=float)
    parser.add_argument("--advance", default=False, help="Advance to next RA")
    args = parser.parse_args()

    if args.initra:
        ra = args.initra
    if args.width:
        width = args.width
    if args.advance:
        advance = args.advance

    for arg in vars(args):
        print(arg, getattr(args, arg))

    while ra < ra_range[1]:
        print("Getting data for RA: ", ra)
        start_time = time.time()
           
        if len(data_obj.tables) == 0:
            parallel_get_data(data_obj, ra, width + overlap)
        
        print("--------------------")
        print("Matching data...")
        
        if advance:
            print("Advancing next RA batch while matching data...")
            tmp_data_obj = TableClass()
            tmp_thread_pool = ThreadPoolExecutor(max_workers=20)
            parallel_get_data_wt_wait(tmp_thread_pool, tmp_data_obj, ra + width, width + overlap)

        matched_data = data_obj.parallel_match_data(previous_df)
        matched_data.to_csv(f'matched_catalog_{ra}.csv', index=False)
        previous_df = matched_data  # Adjust based on expected overlap
        ra += width
        data_obj.reset_tables()

        if advance:
            tmp_thread_pool.shutdown(wait=True)
            data_obj = tmp_data_obj
            tmp_thread_pool = None

        print("--------------------")
        print("Time taken: {:.2f} minutes".format((time.time() - start_time)/60))
        print("--------------------")
    
    
if __name__ == "__main__":
    main()
