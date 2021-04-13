from taigapy import TaigaClient

def get_taiga_id_with_version(taiga_id):
    tc = TaigaClient()
    latest_version = max([int(x['name']) for x in tc.get_dataset_metadata(taiga_id)['versions']])
    taiga_id_with_version = '{}.{}'.format(taiga_id, latest_version)
    return taiga_id_with_version

def get_taiga_ids_list(taiga_id_dict):
    taiga_ids_list = []
    for taiga_id, pairs in taiga_id_dict.items():
        taiga_id_with_version = get_taiga_id_with_version(taiga_id)
        for pair in pairs:
            taiga_ids_list.append({
                'name': pair[0],
                'taiga_id': '{}/{}'.format(taiga_id_with_version, pair[1])
            })
    return taiga_ids_list
