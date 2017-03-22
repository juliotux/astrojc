import os
import appdirs

appname = 'astrojc'

def get_config_dir():
    return appdirs.user_config_dir(appname)
