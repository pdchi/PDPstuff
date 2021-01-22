import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussfit(x, y, double=False, graph=False, quiet=False):
	"""
	Single or double Gaussian fitting routine
	V1.5 for scipy 0.19+

	Parameters:
	x: 1d array of x values to be fit
	y: 1d array of y values to be fit
	double: double or single line fitting, default is single line fitting with double=False
	graph: optional at-a-glance plotting of fitting function
	quiet: remove text response

	Returns:
	6 element 1d array for double
	5 element 1d array for single
    see text output for the definitions of returned values
	"""

	# Set iteration number for curve_fit
	max_fit_repeats = 100000

	# Calculate Initial Gaussian parameters
	peak1 = np.sum(y*x)/np.sum(y)
	FWHM1 = np.sqrt(np.abs(np.sum((x-peak1)**2.*y)/np.sum(y)))
	sig_noise = np.std(y)
	bkg_counts = np.mean(y)
	# Decide if peak is positive or negative w.r.t. continuum
	if np.abs(bkg_counts-np.max(y)) > np.abs(bkg_counts-np.min(y)):
		intensity1 = np.max(y)-bkg_counts
	else:
		intensity1 = np.min(y)-bkg_counts
	peaksep = 2.7
	peakrel = 1.

	## DOUBLE GAUSSIAN FITTING ##
	if double:
		# Define double gaussian curve
		def double_gaussian_fit(x, I, lambda_c, fwhm, R, b):
			return  I*np.exp(-(x-lambda_c)**2.0/(2.0*(fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))))**2.0))+I*np.exp(-(x-lambda_c-R)** 2.0/(2.0*(fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))))**2.0))+b
		# Set parameter bounds and run optimization
		parambounds=([0.,peak1-200,.1,-30,0.],[np.inf,peak1+200,np.inf,30,np.inf])
		popt_double,pcov_double = curve_fit(double_gaussian_fit, x, y, p0=[intensity1, peak1-peaksep, FWHM1*.5, peaksep, bkg_counts], max_nfev=max_fit_repeats,bounds=parambounds)
		# Double Gaussian Chi2
		chi2_double = sum(((y - double_gaussian_fit(x,popt_double[0],popt_double[1],popt_double[2],popt_double[3],popt_double[4]))**2)/ sig_noise**2)  / (len(y) - 5)
		# Equivalent width
		fit = lambda t : popt_double[0]*np.exp(-(t-popt_double[1])**2./(2.*(popt_double[2]/2.355)**2.)) + popt_double[0]*np.exp(-(t-popt_double[1]-popt_double[3])**2./(2.*(popt_double[2]/2.355)**2.))
		eqv = np.sum(-fit(x)) * (x[1] - x[0]) / popt_double[4]
		# Optional text output
		if not quiet:
			print('Fit of double Gaussian, with chi2 of '+str(chi2_double))
			print('Outputs: Peak1 intensity, Peak1 location, FWHM, Relative Sep, Continuum, Equivalent width')
		# Optional graphical output
		if graph:
			plt.clf()
			plt.plot(x,y)
			plt.plot(x,fit(x)+popt_double[4],color='r')
			#plt.axvline(x=popt_double[1]+popt_double[3]+eqv/2.)
			#plt.axvline(x=popt_double[1]+popt_double[3]-eqv/2.)

		return np.append(popt_double,eqv)

	## SINGLE GAUSSIAN FITTING ##
	else:
		# Define single Gaussian curve and optimize
		def single_gaussian_fit(x, I, lambda_c, fwhm, b):
			return  (I*np.exp(-(x-lambda_c)**2.0/(2.0*(fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))))**2.0))) + b
		popt_single,pcov_single = curve_fit(single_gaussian_fit, x, y, p0=[intensity1, peak1, FWHM1, bkg_counts], maxfev=max_fit_repeats)

		# Single Gaussian Chi2
		chi2_single = sum(((y - single_gaussian_fit(x,popt_single[0],popt_single[1],popt_single[2],popt_single[3]))**2) / sig_noise**2)  / (len(y) - 4)

		# Equivalent width
		fit = lambda t : popt_single[0]*np.exp(-(t - popt_single[1])**2./(2. * (popt_single[2]/2.355)**2.0))
		eqv = np.sum(-fit(x)) * (x[1] - x[0]) / popt_single[3]

		# Optional text output
		if not quiet:
			print('Fit of single Gaussian, with chi2 of '+str(chi2_single))
			print('Outputs: Peak intensity, Peak location, FWHM, Continuum, Equivalent Width')

		# Optional graphical output
		if graph:
			plt.clf()
			plt.plot(x,y)
			plt.plot(x,fit(x)+popt_single[3],color='r')
			#plt.axvline(x=popt_single[1]+eqv/2.)
			#plt.axvline(x=popt_single[1]-eqv/2.)

		return np.append(popt_single,eqv)

