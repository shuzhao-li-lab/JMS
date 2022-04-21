import requests
import json
import os

def git_download_from_direcotry(user,repo_name,path_of_dir,local_output_dir):
    '''
    Download files from a github directory

    Input
    =====
    user: user name
    repo_name: the name of repository
    path_of_dir: the path of the github directory (the path after 'tree/master/')
    local_output_dir: the path of the local directory you will to download

    Example
    =====
    git_download_from_direcotry('VirtualMetabolicHuman','AGORA','CurrentVersion/AGORA_1_03/AGORA_1_03_With_Mucins_sbml','../test/input/test_output/')

    '''

    dir_url = os.path.join(f"https://api.github.com/repos/{user}/{repo_name}/contents/",path_of_dir)   
    list_of_attrs =json.loads(requests.get(dir_url).text)
    list_of_dicts = []
    for attrs in list_of_attrs:
        try:
            list_of_dicts.append(
                {'url':attrs['html_url'],
                'name': attrs['name']
                }
            ) # don't use url, which will not work as it go through git API and file size < 1 MB
        except:
            None

    # download the urls
    for dict in list_of_dicts:
        file_path = os.path.join(local_output_dir, dict['name'])
        with open(file_path, 'w') as f:
            url = dict['url']
            r = requests.get(f'{url}?raw=true')
            f.write(r.text)


def git_download_from_file(url, local_output_dir, file_name):
    '''
    Download a github file

    Input
    =====
    url: html_url of the github file
    local_output_dir: the local directory user intended to store the file
    file_name: the intended file name

    Example
    =====
    git_download_from_file("https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.xml","../test/input/test_output/","Human-GEM.xml")

    '''
    file_path = os.path.join(local_output_dir, file_name)
    with open(file_path, 'w') as f:
        r = requests.get(f'{url}?raw=true')
        f.write(r.text)


