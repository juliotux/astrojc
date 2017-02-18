'''
Module to adapt different methods of source detection and background estimation
in standart inputs and outputs.
'''

import numpy as np

class FindAdapterParameters():
    def __init__(self):
        self.detection_fwhm = 3.0
        self.detection_sharplo = 0.0
        self.detection_sharphi = 2.0
        self.detection_roundlo = -1.0
        self.detection_roundhi = 1.0
        self.detection_box_size = 5
        self.background_sigma_clip = 3
        self.background_iters = 1

class FindAdapter():
    def __init__(self, use_photutils=True, use_sep=False,
                 order=['photutils', 'sep']):
        if use_photutils:
            try:
                import photutils
            except:
                use_photutils= False

        if use_sep:
            try:
                import sep
            except:
                use_sep = False

        if not use_photutils and not use_sep:
            raise ImportError('No source finder package available.')

        self.config = FindAdapterParameters()

        self.find_dtype = np.dtype([('x', 'f8'), ('y', 'f8'), ('peak', 'f8'),
                                    ('flux', 'f8'), ('roundness', 'f8'),
                                    ('sharpness', 'f8')])

    def photutils_mmm_background(self, data):
        back = photutils.MMMBackground(photutils.SigmaClip(sigma=self.config.background_sigma_clip,
                                                           iters=self.config.background_iters))
        background = back.calc_background(data)
        bkgrms = photutils.StdBackgroundRMS(photutils.SigmaClip(sigma=self.config.background_sigma_clip,
                                                                iters=self.config.background_iters))
        rms = bkgrms.calc_background_rms(data)

        return background, rms

    def photutils_daofind_source_find(self, data, snr):
        bkg, rms = self.photutils_mmm_background(data)
        thresh = bkg + snr*rms
        daofind = DAOStarFinder(fwhm=self.config.detection_fwhm, threshold=thresh,
                                sharplo=self.config.detection_sharplo,
                                shaprhi=self.config.detection_sharphi,
                                roundlo=self.config.detection_roundlo,
                                roundhi=self.config.detection_roundhi)

        res = daofind(data)
        res_arr = tuple([res[i] for i in ('xcentroid', 'ycentroid', 'peak', 'flux',
                                          'sharpness', 'roundness1')])
        return np.array(list(zip(*res_arr)), dtype=self.find_dtype)

    def photutils_find_peak(self, data):
        bkg, rms = self.photutils_mmm_background(data)
        thresh = bkg + snr*rms
        res = photutils.find_peaks(data, thresh, box_size=self.config.detection_box_size)
        res_arr = tuple([res[i] for i in ('x_peak', 'y_peak', 'peak_value')])
        nan_array = [np.nan]*len(res)

        return np.array(list(zip(*res_arr, nan_array, nan_array, nan_array)),
                        dtype=self.find_dtype)
