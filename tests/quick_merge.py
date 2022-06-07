import json
import sys

if __name__ == "__main__":
    bd = sys.argv[1]
    def read_json(fn):
        return json.load(open(fn))

    out = {"data_object_set": []}

    for dn in ["qa", "assembly", "ReadbasedAnalysis", "annotation", "MAGs"]:
        do = read_json("%s/%s/data_objects.json" % (bd, dn))
        out["data_object_set"].extend(do)
    acts = {
            'annotation': 'metagenome_annotation_activity_set',
            'assembly': 'metagenome_assembly_set',
            'MAGs': 'mags_activity_set',
            'qa': 'read_QC_analysis_activity_set',
            'ReadbasedAnalysis': 'read_based_analysis_activity_set'
            }
    for k, v in acts.items():
        act = read_json("%s/%s/activity.json" % (bd, k))
        act.pop("part_of")
        out[v] = [act]

    print(json.dumps(out))
