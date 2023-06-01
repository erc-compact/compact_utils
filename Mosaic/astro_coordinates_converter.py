import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord

def decimal_to_hours(ra, dec):
    c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    return c.to_string('hmsdms', sep=':', precision=2)

def hours_to_decimal(ra, dec):
    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    return c.ra.deg, c.dec.deg

def equatorial_to_galactic(ra, dec, unit):
    c = SkyCoord(ra, dec, unit=unit)
    l = c.galactic.l.degree
    b = c.galactic.b.degree
    return l, b

def galactic_to_equatorial(l, b, to):
    c = SkyCoord(l, b, frame="galactic", unit=u.deg)
    if to == 'hours':
        return c.to_string('hmsdms', sep=':', precision=2)
    elif to == 'decimal':
        return c.icrs.ra.deg, c.icrs.dec.deg

parser = argparse.ArgumentParser(description="Convert coordinates between equatorial (in hours or decimal degrees) and galactic systems")
parser.add_argument('-ra', type=float, help='Right ascension in decimal degrees or hours')
parser.add_argument('-dec', type=float, help='Declination in decimal degrees')
parser.add_argument('-l', type=float, help='Galactic longitude in decimal degrees')
parser.add_argument('-b', type=float, help='Galactic latitude in decimal degrees')
parser.add_argument('--from', dest='from_', type=str, choices=['hours', 'decimal', 'galactic'], help='The original coordinate system: "hours", "decimal", or "galactic"')
parser.add_argument('--to', type=str, choices=['hours', 'decimal', 'galactic'], help='The conversion target: "hours", "decimal", or "galactic"')

args = parser.parse_args()

if args.from_ == 'hours' and args.to == 'decimal':
    ra_decimal, dec_decimal = hours_to_decimal(args.ra, args.dec)
    print(f'RA: {ra_decimal}, DEC: {dec_decimal}')

elif args.from_ == 'decimal' and args.to == 'hours':
    hmsdms = decimal_to_hours(args.ra, args.dec)
    print(f'RA/DEC in HMS/DMS: {hmsdms}')

elif args.from_ in ['hours', 'decimal'] and args.to == 'galactic':
    l, b = equatorial_to_galactic(args.ra, args.dec, unit=(u.hourangle if args.from_ == 'hours' else u.deg))
    print(f'Galactic l: {l}, b: {b}')

elif args.from_ == 'galactic' and args.to in ['hours', 'decimal']:
    equatorial = galactic_to_equatorial(args.l, args.b, to=args.to)
    print(f'Equatorial coordinates: {equatorial}')
