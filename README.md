# Compton Scattering Experiment

Uses a set of Quantum MCA ASCII data to solve for the mass/energy of an electron and the differential scatting cross section of Cs-137.

Data must be saved in the local directory like so:

```
20deg/20degTI.ASC

20deg/20degTO.ASC

...

100deg/100degTI.ASC

100deg/100degTO.ASC
```
Must be run with python2.7.

## Output

Produces graphs in PNG and SVG format using `pyplot`.

This is the meat of the python script. Here is an A4 letter sized set of plots of the counts over a series of angular changes. This plot is the most intelligent - the peaks are automatically found and a gaussian fit is automatically applied and integrated to find the net counts. The integration regions are shaded in gray and the gaussian fit is composited in smooth blue. The clever observer will note that the trends in the positioning of the peaks as the angle changes roughly corresponds to the data points seen in the Klein Nishina/Thomson scattering differential cross section plot.

![All plots of dataset](https://gitlab.com/charles-m-knox/compton-scattering-experiment/-/raw/master/compton/all_subtracted_plots.png)

Here is an example of the classical/quantum/experimental differential cross sections:

![Klein Nishina formula and Thomson Scattering](https://gitlab.com/charles-m-knox/compton-scattering-experiment/-/raw/master/compton/differential_cross_section_plot.png)

Here we solve for the mass of the electron:

![Mass of electron from experimental data](https://gitlab.com/charles-m-knox/compton-scattering-experiment/-/raw/master/compton/cos_theta_angles.png)
