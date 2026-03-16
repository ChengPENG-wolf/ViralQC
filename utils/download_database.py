import subprocess
import gdown
import os


def download_database(db_path, name):

    print(f'Downloading database to {db_path}')
    if not os.path.exists(db_path):
        os.makedirs(db_path)

    if name == 'refseq':
        id = "1WIcxTJfAFP6BbWYrATsD-tnZG5kzsxzc"
    elif name == 'meta':
        id = '1Z73Y6rJwqIkWjePFxHvmuR6ujN7TNmyp'
    gdown.download(id=id, output=f'{db_path}/viralqc_db.tar.gz', quiet=False)

    cmd = f'tar -xzvf {db_path}/viralqc_db.tar.gz -C {db_path}'
    subprocess.run(cmd, shell=True)
    os.remove(os.path.join(db_path, 'viralqc_db.tar.gz'))

