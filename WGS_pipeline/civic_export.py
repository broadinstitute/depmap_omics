from cravat import constants
from pyliftover import LiftOver
import requests
import json
import time
import csv

def with_retries(attempts, expected_ex, callback):
    for i in range(attempts):
        try:
            result = callback()
        except Exception as ex:
            if isinstance(ex, expected_ex):
                print(f"Got {ex}, retrying in 5 secs")
                time.sleep(5)
                continue
            else:
                raise
        return result
    raise Exception("Too many retry attempts")

def _fetch_civic_table():
    lifter = LiftOver(constants.liftover_chain_paths['hg19'])
    page_url = 'https://civicdb.org/api/variants?count=500&page=1'

    rows = []
    print("Reading from civic...", end="")
    while page_url is not None:
        r = with_retries(20, requests.exceptions.ConnectionError, 
            lambda: requests.get(page_url, timeout=5))
        print(".", end="", flush=True)
        d = json.loads(r.text)
        records = d['records']
        page_url = d['_meta']['links']['next']
        for variant in records:
            row = {}
            if variant['coordinates']['start'] is None or variant['coordinates']['chromosome'] is None:
                continue
            
            row['chromosome_37'] = chrom_37 = variant['coordinates']['chromosome']
            row['start_37'] = pos_37 = int(variant['coordinates']['start'])

            new_coords = lifter.convert_coordinate("chr" + chrom_37, pos_37)
            if len(new_coords) > 0:
                row['chrom'] = new_coords[0][0]#.replace('chr','')
                row['pos'] = new_coords[0][1]
            else:
                continue

            row['ref'] = variant['coordinates']['reference_bases']
            row['alt'] = variant['coordinates']['variant_bases']

            row["description"] = variant.get('description','')
            row["civic_actionability_score"] = variant['civic_actionability_score']
            row["civic_id"] = variant['id']

            rows.append(row)

            # if len(rows) > 10:
            #     break
    print("done", flush=True)
    return rows
        
if __name__ == '__main__':
    rows = _fetch_civic_table()
    with open("civic.csv", "wt") as fd:
        w = csv.DictWriter(fd, rows[0].keys())
        w.writeheader()
        for row in rows:
            w.writerow(row)
