
from astropy.coordinates import Angle


def convert_ra(ra):
    if ' ' in str(ra).strip() or ':' in str(ra).strip():
        ra = str(ra).strip()
        ra = ra.replace(' ', '').replace(':', '')
    else:
        a = Angle(ra, unit='degree').hms
        ra = '{:02.0f}{:02.0f}{:05.2f}'.format(*a)
    return ra


def convert_dec(dec):
    if ' ' in str(dec).strip() or ':' in str(dec).strip():
        dec = str(dec).strip()
        dec = dec.replace(' ', '').replace(':', '')
        if not '-' in dec and not '+' in dec:
            dec = '+' + dec
    else:
        a = Angle(dec, unit='degree').dms
        dec = '{:+03.0f}{:02.0f}{:05.2f}'.format(a[0], abs(a[1]), abs(a[2]))
    return dec


def get_name_from_coordinates(ra, dec):

    if hasattr(ra, '__iter__'):
        ra_ = [convert_ra(x) for x in ra]
        dec_ = [convert_dec(x) for x in dec]
        name = [f'J{r}{d}' for r, d in zip(ra_, dec_)]
    else:
        ra_, dec_ = convert_ra(ra), convert_dec(dec)
        name = f'J{ra_}{dec_}'

    return name