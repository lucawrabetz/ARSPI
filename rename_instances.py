import os
import shutil
import argparse

def copy_directory(A, B):
    # Get the parent directory of A
    parent_directory = os.path.dirname(A)
    
    # Form the absolute paths for directories A and B
    path_A = os.path.abspath(A)
    path_B = os.path.abspath(B)
    
    # Create directory B if it doesn't exist
    os.makedirs(path_B, exist_ok=True)
    
    # Iterate through files in directory A
    for filename in os.listdir(path_A):
        # Check if filename starts with A-
        if filename.startswith(os.path.basename(A) + "-"):
            # Form new filename for directory B
            new_filename = filename.replace(os.path.basename(A), os.path.basename(B), 1)
            
            # Copy the file from A to B
            shutil.copy2(os.path.join(path_A, filename), os.path.join(path_B, new_filename))

def main():
    parser = argparse.ArgumentParser(
        description='Rename an instance set and copy to new or existing directory.')
    parser.add_argument('existing', metavar='X', type=str, nargs=1,
                        help='existing dataframe to add to')
    parser.add_argument('new', metavar='N', type=str, nargs=1,
                        help='existing dataframe to add to')
    args = parser.parse_args()
    old_path = args.existing[0]
    new_path = args.new[0]
    copy_directory(old_path, new_path)

if __name__=='__main__':
    main()
