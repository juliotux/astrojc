import os
import appdirs

appname = 'astrojc'


def get_config_dir():
    return appdirs.user_config_dir(appname)


def get_config_file(filename):
    return os.path.join(get_config_dir(), filename)
