from astropy import units as u

spec_dispersion = u.angstrom

def to_unit(data, unit):
    '''
    Convert data to a unit, avoiding problems with non-quantity data.
    '''
    if hasattr(data, 'unit'):
        return data.to(unit))
    else:
        return data * u.Unit(unit)
