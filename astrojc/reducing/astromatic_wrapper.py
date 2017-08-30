'''
Helpers to run Astromatic software inside python, with results loaded in python.

Partially based on E. Bertin's astromatic_wrapper package.
'''

import os
import subprocess
import traceback
import shutil
import tempfile
from datetime import datetime
import copy
from collections import OrderedDict
import ast
import re

import numpy as np
from six import string_types

from astropy.io import fits, ascii
from astropy.io.votable import parse as parse_vo
from astropy import log as logger
from astropy.table import Table

codes = ['sex', 'psfex', 'swarp', 'scamp', 'skymaker', 'stuff', 'stiff']

class AstromaticError(Exception):
    '''Error for Astromatic realated exceptions.'''

def read_config(f):
    '''Read a astromatic config string (or file.read()).'''
    configs = OrderedDict()
    lines = f.split('\n')
    for i in lines:
        line = i.split()
        #remove comment and blank lines
        if len(line) > 0 and line[0][0] != '#':
            last = len(line)
            for j in range(last):
                #find the last field before a comment
                if line[j][0] == '#':
                    last = j
                    break

            value = ' '.join(line[1:last]).split(',')
            for i in range(len(value)):
                value[i] = value[i].strip(' ')
                try:
                    f = ast.literal_eval(value[i])
                    if isinstance(f, (int, float)):
                        value[i] = f
                except:
                    pass
            if len(value) == 1:
                value = value[0]
            configs[line[0]] = value
    return configs

def write_config(configs):
    '''Write a config dict to a config string in astromatic format.'''
    #find the biggest field name
    big = 0
    for i in configs.keys():
        l = len(str(i))
        if l > big:
            big = l

    config_str = ""
    for i,v in configs.items():
        config_str = config_str + i.ljust(big+1) + '{}'.format(v) + '\n'

    return(config_str)

def find_exe(name, command):
        if command is None:
            command = shutil.which(name)
        try:
            with open(os.devnull, "w") as f:
                subprocess.call(command, stdout=f, stderr=f)
            return command
        except:
            raise

def check_image(image):
    '''Checks if the value passed is a valid file name, an ImageHDU or HDUList.
    If ImageHDU or HDUList is passed, the function will save it as a tmpfile and
    return the name. If a valid file name, the value will be returned. Except,
    None will be returned or an error will be raised.
    '''
    if isinstance(image, string_types):
        if os.path.isfile(image):
            return image
    elif isinstance(image, (fits.PrimaryHDU, fits.ImageHDU)):
        hdulist = fits.HDUList([image])
        tmp = tempfile.mkstemp(prefix='astrojc_', suffix='.fits')
        hdulist.writeto(tmp)
        return tmp
    elif isinstance(image, HDUList):
        tmp = tempfile.mkstemp(prefix='astrojc_', suffix='.fits')
        image.writeto(tmp)
        return tmp
    else:
        raise ValueError('Unrecognized image type.')

class AstromaticProgram():
    '''Base class for run any astromatic program.'''
    def __init__(self, temp_path=None, config={}, config_file=None,
                 store_output=False, store_result=False, **kwargs):
        self.code = None
        self.cmd = None
        self.tmp = temp_path
        self.config = config
        self.config_file = config_file
        self.store_output = store_output
        self.store_result = store_result

    @property
    def version(self):
        return self.get_version()

    @property
    def full_config_list(self):
        '''Returns a full list of the valid configs of SExtractor.'''
        config_str = subprocess.check_output([self._command, '-d']).decode('UTF-8')
        return read_config(config_str).keys()

    def print_config(self):
        print(write_config(self.config))

    def get_cmd(self):
        '''Check if the code is running properly and gen the cmd base list.

        Returns
        -------
        cmd : list
            The list containing the base of the command of the code, like
            ['/usr/bin/sex']
        '''
        if self.cmd is None:
            if self.code not in codes:
                raise AstromaticError('Not valid code.')
            cmd = [find_exe(self.code)]
        else:
            cmd = [self.cmd]

        return cmd

    def get_version(self):
        '''
        Get the version of the currently loaded astromatic code.
        '''
        cmd = self.get_cmd()
        cmd += ['-v']
        try:
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
        except:
            raise AstromaticError('Unable to run `{0}`. Please check that it is'
                                  ' installed correctly'.format(cmd))
        for line in p.stdout.readlines():
            _match = re.search("[Vv]ersion ([0-9\.])+",
                               line.decode(encoding='UTF-8'))
            if _match:
                return str(_match.group()[8:])

    def _build_cmd(self, filenames, **kwargs):
        '''Build a command to run an astromatic code.

        Parameters
        ----------
        filenames: str or list
            Name of a file or list of filenames to run in the command line
            statement
        **kwargs: keyword arguments
            The following are optional keyword arguments that may be used:
            config, config_file
        '''
        # If a single catalog is passed, convert to an array
        if not isinstance(filenames, list):
            filenames = [filenames]

        logger.debug("kwargs used to build command:\n{}".format(kwargs))

        cmd = self.get_cmd()
        cmd += filenames

        # Set a config file to run
        if kwargs.get('config_file', None) is not None:
            cmd += ['-c', kwargs.get('config_file', None)]

        # Put the local configurations in the command line
        local_conf = copy.copy(self.config)
        config = kwargs.get('config', {})
        local_conf.update(config)
        for conf in local_conf.keys():
            if isinstance(local_conf[conf], bool):
                if local_conf[conf]:
                    val = 'Y'
                else:
                    val = 'N'
            else:
                val = local_conf[conf]
            cmd += ['-{}'.format(conf).upper(), str(val)]

        logger.debug('Generated command:\n{}'.format(' '.join(cmd)))

        return cmd, kwargs

    def _run_cmd(self, this_cmd, store_output=False, xml_name=None,
                 raise_error=True, frame=None):
        '''Execute the generated command and handle the outputs.

        Returns
        -------
        result : dict
            Result of the astromatic code execution. This will minimally contain
            a `status` key, that indicates `success` or `error`.
            Additional keys:
            - error_msg : str
                If there is an error and the user is storing the output or
                exporting XML metadata, `error_msg` will contain the error
                message generated by the code
            - output : str
                If `store_output==True` the output of the program execution is
                stored in the `output` value.
            - warnings: str
                If the WRITE_XML parameter is `True` then a table of warnings
                detected in the code is returned
        '''

        result = {'status':'success'}

        # Run code
        logger.info('cmd:\n{0}\n'.format(this_cmd))
        if store_output:
            p = subprocess.Popen(this_cmd, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            p.wait()
            output = p.stdout.read().decode('UTF-8').split('\n')
            result['output'] = output
            for line in output:
                if 'error' in line.lower():
                    result['status'] = 'error'
                    result['error_msg'] = line
                    break
        else:
            status = subprocess.call(this_cmd, shell=True)
            if status > 0:
                result['status'] = 'error'
                if xml_name is not None:
                    votable = parse_vo(xml_name)
                    for param in votable.resources[0].resources[0].params:
                        if param.name == 'Error_Msg':
                            result['error_msg'] = param.value

        # Log any warnings generated by the astromatic code
        if xml_name is not None:
            # SExtractor logs have a '-1' added to the filename and also
            # stream the output catalog to the votable. Since the output may
            # be a FITS_LDAC file, astropy does not rad this properly and it
            # causes the read to crash. This code removes the link to the
            # FITS_LDAC file
            if frame is not None:
                xml_name = xml_name.replace('.xml','-{0}.xml'.format(frame))
            if self.code == 'sex':
                f = open(xml_name, 'r')
                all_lines = f.readlines()
                f.close()
                f = open(xml_name, 'w')
                for line in all_lines:
                    if '<fits' not in line.lower():
                        f.write(line)
                    elif '</DATA>' in line.upper():
                        f.write('</DATA>\n')
                f.close()

            # Sometimes the xml file does not fit the VOTABLE standard,
            # so we mask the invalid parameters
            votable = parse_vo(xml_name, invalid='mask', pedantic=False)
            result['warnings'] = Table.read(votable, table_id='Warnings',
                                            format='votable')
            result['warnings'] = result['warnings'].filled(0)
            result['warnings'].meta['filename'] = xml_name

        if result['status'] == 'error' and raise_error:
            error_msg = "Error in '{0}' execution".format(self.code)
            if 'error_msg' in result:
                error_msg += ': {0}'.format(result['error_msg'])
            raise AstromaticError(error_msg)
        return result

    def run(self, filenames, config={}, config_file=None, **kwargs):
        '''Run the code for some files.'''
        local_conf = copy.copy(self.config)
        local_conf.update(config)

        if config_file is None:
            config_file = self.config_file

        this_cmd, kwargs = self._build_cmd(filenames, config=local_conf,
                                           config_file=config_file,
                                           **kwargs)

        if ('WRITE_XML' in local_conf.keys() and
            'XML_NAME' in local_conf.keys() and
            local_conf['WRITE_XML'] == 'Y'):
            xml_name = local_conf['XML_NAME']
        else:
            xml_name = None

        return self._run_cmd(this_cmd, self.store_output, xml_name)


class SExtractor(AstromaticProgram):
    def __init__(self, cmd_path=None, work_dir=None, config={},
                 config_file=None, store_output=False, store_result=True,
                 **kwargs):
        self.params = kwargs.pop('params', [])
        AstromaticProgram.__init__(self, config=config, config_file=config_file,
                                   store_output=store_output,
                                   store_result=store_result, **kwargs)
        self.code = 'sex'
        if cmd_path is not None:
            self.cmd = cmd_path
        else:
            self.cmd = 'sex'

        if work_dir is None:
            self.work_dir = tempfile.mkdtemp(prefix='astrojc_sex_')
        else:
            try:
                mkdir_p(work_dir)
                self.work_dir = work_dir
            except:
                self.work_dir = tempfile.mkdtemp(prefix='astrojc_sex_')
                logger.warn('Could not create the given work_dir, using tmp.')

    @property
    def full_params_list(self):
        '''Return the full list of SExtractor parameters.'''
        full_params = subprocess.check_output([self.cmd, '-dp']).decode('UTF-8')
        return list(map(lambda s: s[1:-1],
                        re.compile("#\w*\s").findall(full_params)))

    def _check_params(self, params):
        '''Check if the given parameters are in the full_params_list.'''
        full_params_list = self.full_params_list
        for i in params:
            if i not in full_params_list:
                return False
        return True

    def _read_result(self, cat_name, cat_type):
        '''Read the sextractor result.'''

    def run(self, image, flag_img=None, weight_img=None, params=None,
            cat_type='ASCII_HEAD', cat_name=None, config={}, config_file=None,
            **kwargs):
        '''Run sextractor for a image.

        Parameters
        ----------
        image : str or astropy.io.fits.ImageHDU
            Image name or a fits HDU containing the image to process.
        flag_img : str or astropy.io.fits.ImageHDU
            Image name or a fits HDU containing the quality flags.
        weight_img : str or astropy.io.fits.ImageHDU
            Image name or a fits HDU containing the weight map.
        params : list_like
            List of the parameters to be outputed.
        cat_type : str
            The format of output catalog, to be placed in `CATALOG_TYPE` in the
            configure. Default: `ASCII_HEAD`. Override any previous definition.
        cat_name : str
            The name of the output catalog, to be placed in `CATALOG_NAME` key
            in the configure. If None, no output catalog will be saved.
        config : dict_like
            Dictionary of configs to be used in the run. Override previous
            configure items.
        config_file : str
            Custom name of a config file to be passed to `sex`. If None,
            the default config file will be tryed, except, no file will be
            passed and `sex` will use defaults.
        kwargs
            Additional kwargs to be passed to AstromaticProgram.run()
        '''
        local_conf = copy.copy(self.config)
        local_conf.update(config)

        if config_file is None:
            config_file = self.config_file

        if params is not None:
            if 'PARAMETERS_NAME' in local_conf.keys():
                logger.warn('Multiple parameter files specified, using local '
                            '`params` list.'.format(params))
            logger.debug('Using params list:\n{}'.format(params))
            if not self._check_params(params):
                logger.warn('Unrecognized param in list.')
            param_name = os.path.join(self.work_dir, 'sex.param')
            f = open(param_name, 'w')
            f.write('\n'.join(params))
            f.close()
            local_conf['PARAMETERS_NAME'] = param_name
        elif 'PARAMETERS_NAME' not in local_conf.keys():
            raise AstromaticError('To run SExtractor yo must either supply a '
                                  '`params` list of parameters or a config '
                                  'keyword `PARAMETERS_NAME` with the .param '
                                  'file.')

        image_name = check_image(image)
        if flag_img is not None:
            self.config['FLAG_IMAGE'] = check_image(flag_img)
        if weight_img is not None:
            self.config['WEIGHT_IMAGE'] = check_image(weight_img)

        cat_type = local_conf.get('CATALOG_TYPE', cat_type)
        cat_name = local_conf.get('CATALOG_NAME', cat_name)
        if cat_name is not None:
            self.config['CATALOG_NAME'] = cat_name
        else:
            if cat_type in ['FITS_1.0', 'FITS_LDAC']:
                ext = 'fits'
            else:
                ext = 'cat'
            self.config['CATALOG_NAME'] = os.path.join(self.work_dir,
                                                       'result.{}'.format(ext))
        self.config['CATALOG_TYPE'] = cat_type

        stat = AstromaticProgram.run(self, image, config=local_conf,
                                     config_file=config_file, **kwargs)

        if stat['status'] == 'success' and self.store_result:
            result = self._read_result(cat_name, cat_type)
        else:
            result = None

        if self.store_result:
            stat['data'] = result

        return stat

def sex(image, params=[], configs={}, sex_command=None, work_dir=None):
    '''Run sextractor for a image.

    Returns a table with default+given parameters list and default+given configs
    dictionary.
    The sextractor command can be set by the sex_command variable.

    Parameters:
        image : array_like or fits.HDUList or fits.ImageHDU
            The image to run in sextractor.
        params : list_like
            The list of the parameters to the sextractor output.
        configs : dict_like
            Custom configures to be updated in the default config dict.
        sex_command : string
            The `sex` command to run. If None, we will locate it with `which`.
            Default: None
        work_dir : string
            If set, the files will be saved in this folder. Except, a tmp dir
            will be created.
    '''
    #Now, just for testing
    sex_default_params = ['XWIN_IMAGE', 'YWIN_IMAGE', 'AWIN_IMAGE', 'BWIN_IMAGE',
                          'THETAWIN_IMAGE', 'BACKGROUND', 'FLUX_AUTO']
    s = SExtractor(store_output=True, store_result=True,
                   config={'FILTER' : 'N'})
    print(s.run('/home/julio/17fev21/AGCar_A3_L0_F2.fits', cat_type='FITS_LDAC',
                params = sex_default_params))
    #print(s.default_config)
    #print(s.default_params)




sex(None)
