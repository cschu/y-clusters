import json
import re
import sys

import pymongo

import pprint as pp


with open(sys.argv[1]) as _in:
    conn_data = json.load(_in)

records = []


if sys.argv[2] == "spire":
    conn_data = conn_data['mags']
    client = pymongo.MongoClient(f"mongodb://{conn_data['user']}:{conn_data['password']}@{conn_data['host']}:{conn_data['port']}")

    db = client['mags']

    # we want all binned MAGs (=bins) with existing spire_v1_cluster (=specI cluster)
    for record in db.bins.find(
        {
            '$and': [
                {'spire_v1_cluster': {'$exists': True}},
                {
                    '$or': [
                        {'spire_v1_cluster.assignment_method': 'pg_v3_mapped_ANI_95'},
                        {'spire_v1_cluster.assignment_method': 'pg_v3_mapped_marker_gene'},
                    ]
                }
            ]
        }
    ):
        # print(record.get('spire_v1_cluster', {}).get('name'), record.get('bin_id'), record.get('bin_path'), sep='\t')
        print(record.get('bin_id'), re.sub(r'specI_v4_0*([0-9]+)', lambda x:x.group(1) , record.get('spire_v1_cluster', {}).get('name', '')), sep='\t')

elif sys.argv[2] == "pg3":
    conn_data = conn_data['progenomes']
    client = pymongo.MongoClient(f"mongodb://{conn_data['user']}:{conn_data['password']}@{conn_data['host']}:{conn_data['port']}")

    db = client['progenomes']

    # we only want isolate genomes with assigned fr13_cluster (=specI cluster)
    for record in db.samples.find({'fr13_cluster': {'$exists': True}}):
        print(record.get('fr13_cluster'), record.get('sample_id'), record.get('analysis_path'), sep="\t")
