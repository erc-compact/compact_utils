
import numpy as np
import pandas as pd
import sqlalchemy
from sqlalchemy import text
import configparser
import uuid
import base64
import os

database_ini = os.environ.get('HOME') + '/.compactdb_secrets.ini'

config = configparser.ConfigParser()
config.read(database_ini)

DB_USER = config['compactdb_ro']['DB_USER']
DB_PASSWORD = config['compactdb_ro']['DB_PASSWORD']
DB_HOST = config['compactdb_ro']['DB_HOST']
DB_PORT = config['compactdb_ro']['DB_PORT']
DB_NAME = config['compactdb_ro']['DB_NAME']

DATABASE_URI = (f"mysql+pymysql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}")
engine = sqlalchemy.create_engine(DATABASE_URI)


class UUIDUtility:
    @staticmethod
    def generate_uuid_string():
        return str(uuid.uuid4())

    @staticmethod
    def generate_binary_uuid():
        return uuid.uuid4().bytes

    @staticmethod
    def convert_uuid_string_to_binary(uuid_string):
        return uuid.UUID(uuid_string).bytes

    @staticmethod
    def convert_binary_uuid_to_base64(binary_uuid):
        return base64.b64encode(binary_uuid).decode('ascii')

    @staticmethod
    def convert_binary_uuid_to_string(binary_uuid):
        # Ensure binary_uuid is of type bytes
        if isinstance(binary_uuid, bytearray):
            binary_uuid = bytes(binary_uuid)
        return str(uuid.UUID(bytes=binary_uuid))
    
    @staticmethod
    def generate_uuid_list(n):
        # Each UUID is 16 bytes.
        # os.urandom(...) is cryptographically secure in Python 2.7.
        random_bytes = os.urandom(16 * n)
        return [
            str(uuid.UUID(bytes=random_bytes[i*16:(i+1)*16])) 
            for i in range(n)
        ]
    
shortlist_query = """
SELECT 
    unique(HEX(f.id)) AS fold_candidate_id
FROM fold_candidate f
INNER JOIN candidate_tracker c
    ON f.id = c.fold_candidate_id
INNER JOIN  
    data_product dp ON dp.id = f.dp_id
WHERE
    dp.`created_by_run_name` IN ("CRAZY_GUS_4", "CRAZY_GUS_5", "CRAZY_GUS_6", "CRAZY_GUS_7", "CRAZY_GUS_8","CRAZY_GUS_9","CRAZY_GUS_10","CRAZY_GUS_11") and
	( f.fold_snr > 7.5 OR
      (c.candidate_filter_id = 1 AND c.value > 0.1)
   OR
      (c.candidate_filter_id = 2 AND c.value > 0.1)
   OR
      (c.candidate_filter_id = 3 AND c.value > 0.1)
   OR
      (c.candidate_filter_id = 4 AND c.value > 0.1) )
"""

shortlist_df = pd.read_sql_query(shortlist_query, engine)

print(f"Shortlisted data frame has {len(shortlist_df.index)} number of rows")


main_query = f"""
    SELECT 
        po.id as pointing_id, 
        b.id as beam_id, 
        b.name as beam_name, 
        t.target_name as source_name,
        b.ra_str as ra, 
        b.dec_str as "dec", 
        '0' as gl, 
        '0' as gb, 
        '0' as mjd_start, 
        po.utc_start as utc_start, 
        1/s.spin_period as f0_usr, 
        1/f.spin_period as f0_opt, 
        '0' as f0_opt_err, 
        -1 * s.pdot / (s.spin_period  * s.spin_period) as f1_user, 
        -1 * f.pdot / (f.spin_period  * f.spin_period) as f1_opt, 
        '0' as f1_opt_err, 
        s.pdot * 3e8 / s.spin_period as acc_user, 
        f.pdot * 3e8 / f.spin_period as acc_opt, 
        '0' as acc_opt_err, 
        s.dm as dm_user, 
        f.dm as dm_opt, 
        '0' as dm_opt_err, 
        s.snr as sn_fft, 
        f.fold_snr as sn_fold, 
        s.segment_pepoch as pepoch, 
        '0' as maxdm_ymw16, 
        '0' as dist_ymw16,
        HEX(f.id) AS fold_candidate_id, 
        f.id as fold_candidate_id_bin,
        CONCAT(dp.filepath, CONCAT(LEFT(dp.filename, LENGTH(dp.filename) - 3), '.png')) AS png_path,
        "null" as filterbank_path, 
        "metafiles/3HM.meta" AS metafile_path,
        "null" as candidate_tarball_path


FROM fold_candidate f
INNER JOIN search_candidate s ON s.id = f.search_candidate_id
INNER JOIN data_product dp ON f.dp_id = dp.id
INNER JOIN processing p ON dp.created_by_task_id = p.task_id AND dp.created_by_run_name = p.run_name
INNER JOIN beam b ON b.id = dp.beam_id
INNER JOIN pointing po ON po.id = b.pointing_id
INNER JOIN target t ON t.id = po.target_id
WHERE HEX(f.id) IN ({", ".join([f"'{id_}'" for id_ in shortlist_df['fold_candidate_id'].tolist()])})
"""
main_df = pd.read_sql_query(main_query, engine)

print(f"Main DF has {len(main_df.index)} number of rows")


# 2) Extract IDs as a list
id_list = main_df["fold_candidate_id"].tolist()

# 3) Build placeholders & param dictionary for the IN clause
placeholders = []
params = {}
for i, val in enumerate(id_list):
    key = f"id_{i}"
    placeholders.append(f":{key}")
    params[key] = val

# 4) Example pics_palfa query
pics_palfa_query = text(f"""
SELECT 
    c.value AS pics_palfa,
    HEX(f.id) AS fold_candidate_id
FROM fold_candidate f
INNER JOIN candidate_tracker c 
    ON f.id = c.fold_candidate_id
WHERE c.candidate_filter_id = 1
  AND HEX(f.id) IN ({", ".join(placeholders)})
""")

pics_trapum_ter5_query = text(f"""
SELECT
    c.value AS pics_trapum_ter5,
    HEX(f.id) AS fold_candidate_id
FROM fold_candidate f
INNER JOIN candidate_tracker c 
    ON f.id = c.fold_candidate_id
WHERE c.candidate_filter_id = 2
  AND HEX(f.id) IN ({", ".join(placeholders)})
""")

pics_palfa_meerkat_l_sband_best_fscore = text(f"""
SELECT
    c.value AS pics_palfa_meerkat_l_sband_best_fscore,
    HEX(f.id) AS fold_candidate_id
FROM fold_candidate f
INNER JOIN candidate_tracker c 
    ON f.id = c.fold_candidate_id
WHERE c.candidate_filter_id = 3
  AND HEX(f.id) IN ({", ".join(placeholders)})
""")


pics_meerkat_l_sband_combined_best_recall = text(f"""
SELECT
    c.value AS pics_meerkat_l_sband_combined_best_recall,
    HEX(f.id) AS fold_candidate_id
FROM fold_candidate f
INNER JOIN candidate_tracker c 
    ON f.id = c.fold_candidate_id
WHERE c.candidate_filter_id = 4
  AND HEX(f.id) IN ({", ".join(placeholders)})
""")

pics_palfa_df = pd.read_sql_query(pics_palfa_query, engine, params=params)
pics_trapum_ter5_df = pd.read_sql_query(pics_trapum_ter5_query, engine, params=params)
pics_palfa_meerkat_l_sband_best_fscore_df = pd.read_sql_query(pics_palfa_meerkat_l_sband_best_fscore, engine, params=params)
pics_meerkat_l_sband_combined_best_recall_df = pd.read_sql_query(pics_meerkat_l_sband_combined_best_recall, engine, params=params)

                                  
merged_df = main_df.merge(pics_palfa_df, on="fold_candidate_id", how="left")
merged_df = merged_df.merge(pics_trapum_ter5_df, on="fold_candidate_id", how="left")
merged_df = merged_df.merge(pics_palfa_meerkat_l_sband_best_fscore_df, on="fold_candidate_id", how="left")
merged_df = merged_df.merge(pics_meerkat_l_sband_combined_best_recall_df, on="fold_candidate_id", how="left")

#convert fold_candidate_id_bin to string using uuid
merged_df['fold_candidate_id_bin'] = merged_df['fold_candidate_id_bin'].apply(UUIDUtility.convert_binary_uuid_to_string)

#remove fold_candidate_id column
merged_df = merged_df.drop(columns=['fold_candidate_id'])

#rename png_path as original_png_path
merged_df.rename(columns={'png_path':'original_png_path'}, inplace=True)

#write out png_path as a new column which is plots/{basename(png_path)}
merged_df['png_path'] = 'plots/' + merged_df['original_png_path'].str.split('/').str[-1]

print(f"merged df has {len(merged_df.index)} rows")



#save to csv
merged_df.to_csv('merged_df_round2.csv', index=False)




