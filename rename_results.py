import os
import pandas as pd
import argparse

def change_setname(new_base):
    return new_base

def change_instancename(instance_name, new_base):
    graph_params = instance_name.split('-')[1]
    cost_params = instance_name.split('-')[2]
    return new_base + '-' + graph_params + '-' + cost_params

def main():
    parser = argparse.ArgumentParser(
        description='Rename all instances and setnames in a results file.')
    parser.add_argument('input', metavar='I', type=str, nargs=1,
                        help='existing results file path')
    parser.add_argument('new', metavar='N', type=str, nargs=1,
                        help='new results file path')
    args = parser.parse_args()
    input_file = args.input[0]
    output_file = args.new[0]
    new_base = os.path.basename(output_file).split('.')[0]
    df = pd.read_csv(input_file)
    df['set_name'] = df['set_name'].apply(lambda _: change_setname(new_base))
    df['instance_name'] = df['instance_name'].apply(lambda col : change_instancename(col, new_base))
    df.to_csv(output_file, header=True, index=False)

if __name__=='__main__':
    main()
