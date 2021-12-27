# This is a sample Python script.
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

import matplotlib.pyplot as plt
# %matplotlib inline
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.io import fits
from astropy.table import QTable
from astropy.utils.data import download_file

from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = 10000  # Set the row limit for returned data

#2
# ngc188_center = SkyCoord(12.11*u.deg, 85.26*u.deg)
# ngc188_center

#3
# ngc188_center = SkyCoord(12.11*u.deg, 85.26*u.deg, frame='icrs')
# ngc188_center

#6
ngc188_center = SkyCoord.from_name('NGC 188', frame='icrs')
ngc188_center

job = Gaia.cone_search_async(ngc188_center, radius=0.5*u.deg)
ngc188_table = job.get_results()

# only keep stars brighter than G=19 magnitude
ngc188_table = ngc188_table[ngc188_table['phot_g_mean_mag'] < 19*u.mag]

#15
cols = ['source_id',
 'ra',
 'dec',
 'parallax',
 'parallax_error',
 'pmra',
 'pmdec',
 'radial_velocity',
 'phot_g_mean_mag',
 'phot_bp_mean_mag',
 'phot_rp_mean_mag']

#16
ngc188_table[cols].write('gaia_results.fits', overwrite=True)

#17
print(len(ngc188_table))

#20
ngc188_gaia_coords = SkyCoord(ngc188_table['ra'], ngc188_table['dec'])

#24
ngc188_center_3d = SkyCoord(12.11*u.deg, 85.26*u.deg, distance=1.96*u.kpc)

#25
parallax_snr = ngc188_table['parallax'] / ngc188_table['parallax_error']
ngc188_table_3d = ngc188_table[parallax_snr > 10]
print(len(ngc188_table_3d))

#26
# Distance(parallax=1*u.mas)

#27
gaia_dist = Distance(parallax=ngc188_table_3d['parallax'])

#28
ngc188_coords_3d = SkyCoord(ra=ngc188_table_3d['ra'],
                            dec=ngc188_table_3d['dec'],
                            distance=gaia_dist)

#29
fig, ax = plt.subplots(figsize=(6.5, 5.2),
                       constrained_layout=True)
cs = ax.scatter(ngc188_coords_3d.ra.degree,
                ngc188_coords_3d.dec.degree,
                c=ngc188_coords_3d.distance.kpc,
                s=5, vmin=1.5, vmax=2.5, cmap='twilight')
cb = fig.colorbar(cs)
cb.set_label(f'distance [{u.kpc:latex_inline}]')

ax.set_xlabel('RA [deg]')
ax.set_ylabel('Dec [deg]')

ax.set_title('Gaia DR2 sources near NGC 188', fontsize=18)

# display on Raspberry Pi?

plt.savefig('cool_plot-2.png')

from kivy.app import App
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label
from kivy.uix.image import Image
from kivy.uix.button import Button
from kivy.uix.textinput import TextInput

class DisplayPlot(App):
    def build(self):
        self.window = GridLayout()
        self.window.cols = 1

        self.window.add_widget(Image(source="cool_plot-2.png"))

        return self.window

DisplayPlot().run()