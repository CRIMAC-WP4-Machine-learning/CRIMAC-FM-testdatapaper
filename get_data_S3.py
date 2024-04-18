import os
import boto3
from tqdm import tqdm


""""
Copyright 2021 the Norwegian Computing Center

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
"""

# This script dosnloads the data from crimac-scratch/CRIMAC-FM-testdata to your local
# copy of the $CRIMACSCRATCH.


def boto3download(host, access_key, secret_key, bucketname, s3folder, savefolder):
    s3 = boto3.resource('s3', aws_access_key_id=access_key, aws_secret_access_key=secret_key, region_name='us-east-1',
                        endpoint_url=host, )
    bucket = s3.Bucket(bucketname)

    objects = bucket.objects.filter(Prefix=s3folder)

    pbar = tqdm(objects)
    for my_bucket_object in pbar:
        savefile = my_bucket_object.key.replace("gpfs0-crimac-scratch/", "")
        path = os.path.join(savefolder, savefile)
        dirname = os.path.dirname(path) + os.sep
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        # print(path)
        pbar.set_description(f"Downloading {path}")
        bucket.download_file(my_bucket_object.key, path)


def show_folders(host, access_key, secret_key, bucketname, s3folder):
    """
    Function to show the directory tree in a S2 bucket
    :param host:
    :param access_key:
    :param secret_key:
    :param bucketname:
    :param s3folder:
    :return:
    """
    s3 = boto3.resource('s3', aws_access_key_id=access_key, aws_secret_access_key=secret_key, region_name='us-east-1',
                        endpoint_url=host, )

    bucket = s3.Bucket(bucketname)

    objects = []
    for bucket_object in bucket.objects.filter(Prefix=s3folder):
        objects.append(bucket_object.key)

    tree = {}
    for obj in objects:
        split = obj.split('/')

        cur_tree = tree
        for subdir in split:
            if '.zarr' in subdir:
                cur_tree[subdir] = {}
                break
            if subdir not in cur_tree.keys():
                cur_tree[subdir] = {}

            cur_tree = cur_tree[subdir]

    space = '    '
    branch = '│   '
    tee = '├── '
    last = '└── '

    def print_tree(tree, spaces=''):
        new_space = branch + spaces
        for key in tree.keys():
            print(spaces + tee + key)
            print_tree(tree[key], spaces=new_space)

    print()
    print_tree(tree)
    print()


host = 'https://s3.hi.no'
access_key = "crimac"  # Username
secret_key = "9!%L*h7Q"  # Password
bucketname = 'crimac-scratch'  # s3 bucket

s3folder = 'gpfs0-crimac-scratch/CRIMAC-FM-testdata/'  # all s3  crimac folders at HI start with gpfs0-crimac-scratch/
savefolder = os.path.join(os.getenv('CRIMACSCRATCH'), 'CRIMAC-FM-testdata')

# for s3folder in data:
show_folders(host, access_key, secret_key, bucketname, s3folder)
print("Downloading data from s3 bucket")
boto3download(host, access_key, secret_key, bucketname, s3folder, savefolder)

