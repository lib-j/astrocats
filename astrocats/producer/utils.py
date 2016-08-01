"""
"""
import os
import requests


def event_filename(name):
    return (name.replace('/', '_'))


def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)


def utf8(x):
    return str(x, 'utf-8')


def label_format(label):
    newlabel = label.replace('Angstrom', 'Å')
    newlabel = newlabel.replace('^2', '²')
    return newlabel


def is_valid_link(url):
    response = requests.get(url)
    try:
        response.raise_for_status()
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        return False
    return True


def get_first_value(catalog, name, field):
    return catalog[name][field][0]['value'] if field in catalog[
        name] and catalog[name][field] else ''


def get_first_kind(catalog, name, field):
    return (catalog[name][field][0]['kind'] if field in catalog[name] and
            catalog[name][field] and 'kind' in catalog[name][field][0] else '')


def md5file(fname):
    import hashlib
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
