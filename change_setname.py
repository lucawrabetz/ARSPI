import pandas as pd

INPUT = 'follower_testbed_final.csv'
NEW_SETNAME = 'experiment_bendersenum'
OUTFILE = NEW_SETNAME + '.csv'


def change_setname(set_name):
    return NEW_SETNAME


def change_instancename(instance_name):
    graph_params = instance_name.split('-')[1]
    cost_params = instance_name.split('-')[2]
    return NEW_SETNAME + '-' + graph_params + '-' + cost_params


df = pd.read_csv(INPUT)
df['set_name'] = df['set_name'].apply(change_setname)
df['instance_name'] = df['instance_name'].apply(change_instancename)
df.to_csv(OUTFILE, header=True, index=False)
