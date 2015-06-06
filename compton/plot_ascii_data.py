import matplotlib.pyplot as plt
import argparse
import collections
import matplotlib.mlab as mlab
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import math

def main(args_dict):
	def gauss_function(x, a, x0, sigma):
		return a*np.exp(-(x-x0)**2/(2*sigma**2))

	def wavelength_shift_function(E_gamma, theta_degrees, m_e=.511):
		return 1/(1+(E_gamma/(m_e)*(1-math.cos(math.radians(theta_degrees)))))

	def klein_nishina_function(E_gamma, theta_degrees, alpha_const=1/137.04, r_c_const=2.42631023*10**(-12)/(2*3.14159),m_e=.511):
		wavelength_result = wavelength_shift_function(E_gamma, theta_degrees, m_e)
		return (alpha_const**2)*(r_c_const**2)*(wavelength_result**2)*(wavelength_result+wavelength_result**(-1)-1+(math.cos(math.radians(theta_degrees))**2))/2

	def get_x_y_values_from_data_lines(data_lines):
		#Skip the first 12 lines
		x_values = []
		y_values = []
		for x in range(11, len(data_lines)):
			data_points_in_line = []
			for s in data_lines[x].split():
				if s.isdigit():
					data_points_in_line.append(s)
			if len(data_points_in_line) == 2:
				x_values.append(data_points_in_line[0])
				y_values.append(data_points_in_line[1])
		return x_values, y_values

	def get_calibrations_from_text_lines(input_data_lines):
		calibration_A = float(input_data_lines[6].split("A =")[1].strip())
		calibration_B = float(input_data_lines[7].split("B =")[1].strip())
		calibration_C = float(input_data_lines[8].split("C =")[1].strip())
		return calibration_A,calibration_B,calibration_C

	def get_peak_value_simple(all_input_data, low_energy_threshold=200, high_energy_threshold=625):
		#Picks out the max value for a given set of energy/bins
		#Returns the peak bin number and the energy value it occurs at
		#low_energy_threshold = value to ignore all values below
		peak_bins = 0
		peak_energy = 0
		for datum in all_input_data:
			if int(all_input_data[datum]) > peak_bins and int(datum) > low_energy_threshold and int(datum) < high_energy_threshold:
				peak_bins = all_input_data[datum]
				peak_energy = datum
		return int(peak_bins), int(peak_energy)

	def get_average_values_from_data(all_input_data, sample_size=200, energy_min=200, energy_max=800):
		#Prints average values (trying to find what the average value for ~0 is)
		current_datapoint_index = 0
		current_sample_total = 0
		current_sample_average = 0
		for datum in all_input_data:
			if datum > energy_min and datum < energy_max:
				if current_datapoint_index % sample_size == 0: #done with sample
					current_sample_average = current_sample_total / sample_size
					#print "avg energy over interval: x=("+str(datum-sample_size)+","+str(datum)+") = " + str(current_sample_average)
					#plt.annotate(s=str(int(current_sample_average)), xy=(datum, 300) )
					current_sample_average = 0
					current_sample_total = 0
				else:
					current_sample_total += all_input_data[datum]
				current_datapoint_index += 1

	def get_gaussian_peak_width(all_input_data, peak_energy, sample_size=int(args_dict['peak_width_sample_size']), energy_min=200, energy_max=800, avg_threshold=int(args_dict["averaging_threshold"])):
		#On average, the width of each peak is about 100-150 KeV
		#The negligible data has an average value of about 0-28 => samples below the avg_threshold will be considered 0
		#start from the gaussian peak - then get the average value of 50 samples => is avg below threshold? { No? keep going
		#...................................................................................................{ Yes? return energy
		#Returns the width of the peak!
		current_datapoint_index = 0
		current_sample_total = 0
		current_sample_average = 0
		#Calculate Right Side first
		for datum in all_input_data:
			if datum > peak_energy:
				if current_datapoint_index % sample_size == 0:
					current_sample_average = current_sample_total / sample_size
					peak_width = (datum - peak_energy) * 2
					if current_sample_average > avg_threshold:
						return peak_width
					current_sample_average = 0
					current_sample_total = 0
				else:
					current_sample_total += all_input_data[datum]
				current_datapoint_index += 1
			#else:
			#	#Ignore this datapoint
		#Then do left side

	def get_all_plots():
		plt.figure(figsize=(8.5,22))
		plot_index = 0
		for i in range(20, 110, 10):
			plot_index += 1
			#Target Out Plot
			current_file_name_target_out = str(i) + "deg/"+str(i)+"degTO.ASC"
			input_file = open(current_file_name_target_out, "r")
			input_text_lines = input_file.readlines()
			channel_vals, bin_vals = get_x_y_values_from_data_lines(input_text_lines)
			cal_A, cal_B, cal_C = get_calibrations_from_text_lines(input_text_lines)
			energy_vals = []
			for j in range(0, len(channel_vals)):
				energy_vals.append(float(cal_A + cal_B*float(channel_vals[j]) + (cal_C*float(channel_vals[j])**2)))
			plt.subplot(9,2, plot_index)
			plt.title(current_file_name_target_out)
			#plt.fill(energy_vals, bin_vals,'b')
			plt.plot(energy_vals, bin_vals)

			#Target In Plot
			plot_index += 1
			current_file_name_target_in = str(i) + "deg/"+str(i)+"degTI.ASC"
			input_file = open(current_file_name_target_in, "r")
			input_text_lines = input_file.readlines()
			channel_vals, bin_vals = get_x_y_values_from_data_lines(input_text_lines)
			cal_A, cal_B, cal_C = get_calibrations_from_text_lines(input_text_lines)
			energy_vals = []
			for j in range(0, len(channel_vals)):
				energy_vals.append(float(cal_A + cal_B*float(channel_vals[j]) + (cal_C*float(channel_vals[j])**2)))
			plt.subplot(9,2, plot_index)
			plt.title(current_file_name_target_in)
			#plt.fill(energy_vals, bin_vals,'b')
			plt.plot(energy_vals, bin_vals)
		plt.savefig("all_plots_out.svg")
		plt.savefig("all_plots_out.png")

	def get_subtracted_plots():
		angles_energies = []
		plt.figure(figsize=(8.5,11))
		plot_index = 0
		for i in range(20, 110, 10):
			plot_index += 1
			#Target Out Plot
			current_file_name_target_out = str(i) + "deg/"+str(i)+"degTO.ASC"
			input_file = open(current_file_name_target_out, "r")
			input_text_lines = input_file.readlines()
			cal_A, cal_B, cal_C = get_calibrations_from_text_lines(input_text_lines)
			channels_out, bin_out = get_x_y_values_from_data_lines(input_text_lines)
			energies_out = []
			for j in range(0,len(channels_out)):
				energies_out.append(float(cal_A + cal_B*float(channels_out[j]) + (cal_C*float(channels_out[j])**2)))
			#Target In Plot
			current_file_name_target_in = str(i) + "deg/"+str(i)+"degTI.ASC"
			input_file = open(current_file_name_target_in, "r")
			input_text_lines = input_file.readlines()
			channels_in, bin_in = get_x_y_values_from_data_lines(input_text_lines)
			energies_in = []
			for j in range(0,len(channels_in)):
				energies_in.append(float(cal_A + cal_B*float(channels_in[j]) + (cal_C*float(channels_in[j])**2)))

			#Plot x = energy, y = bin_in - bin_out
			yvals = []
			xvals = []

			all_points = {}
			for x in range(0,len(channels_in)):
				xval = (energies_out[x] + energies_in[x])/2 #Take an average of the energy at each point
				yval = float(bin_in[x]) - float(bin_out[x])
				xvals.append(xval) 
				yvals.append(yval)
				all_points[xval] = yval
			all_points_sorted = collections.OrderedDict(sorted(all_points.items(), key=lambda t: t[0]))

			plt.subplot(9,1, plot_index)
			plt.title(str(i) + "Deg TI - TO")
			plt.xlim(150,800)
			plt.ylim(-200,1000)
			peak_bins, peak_energy = get_peak_value_simple(all_points_sorted)
			#plt.text(str(peak_energy), str(peak_bins), "("+str(peak_energy)+","+str(peak_bins)+")")
			plt.xticks([200,400,600,800])
			plt.yticks([0,200,-200])
			plt.ylabel('Counts')
			plt.xlabel('KeV')
			plt.grid(axis="y", alpha=0.25, linestyle="-")
			#plt.fill(xvals, yvals,'b')
			#get_average_values_from_data(all_points_sorted)
			this_peak_width = str(get_gaussian_peak_width(all_points_sorted, peak_energy))
			#plt.annotate(s="Avg: "+str(float(this_peak_width))+", Width: " + str(peak_width), xy=(datum, 300))
			true_energy_vals = all_points_sorted.keys()
			true_bin_vals = all_points_sorted.values()
			#Get energy vals for plotting normal distribution
			normal_dist_energies = []
			normal_dist_bins = []
			for datum in all_points_sorted:
				if datum > (float(peak_energy) - float(this_peak_width)) and datum < (float(peak_energy) + float(this_peak_width)):
					normal_dist_energies.append(datum)
					normal_dist_bins.append(all_points_sorted[datum])
			#estimate mean and standard deviation
			mean = np.mean(np.array([normal_dist_energies, normal_dist_bins]),dtype=np.float64)
			sigma = np.std(np.array([normal_dist_energies, normal_dist_bins]),dtype=np.float64)
			#do the fit!
			popt, pcov = curve_fit(gauss_function, normal_dist_energies, normal_dist_bins, p0 = [1, mean, sigma])
			#popt[1] is the ACTUAL gaussian peak!! => peak energy!
			#popt[0] is the actual gaussian width
			#Do the last few lines again with the new results
			peak_energy = popt[1]
			peak_half_width = popt[0] / 2
			normal_dist_energies = []
			normal_dist_bins = []
			for datum in all_points_sorted:
				if datum > (float(peak_energy) - float(peak_half_width)) and datum < (float(peak_energy) + float(peak_half_width)):
					normal_dist_energies.append(datum)
					normal_dist_bins.append(all_points_sorted[datum])
			plt.axvspan(float(peak_energy) - float(peak_half_width), float(peak_energy)+float(peak_half_width), color='gray', alpha=0.25)

			#plot the fit results
			mean_neat = int(mean)
			sigma_neat = int(sigma)
			y_top = 1000
			#Integrate: total up the count over our peak width
			total_counts_from_gaussian = []
			#get baseline offset from horizontal axis for gaussian to reduce bad data
			#vertical_offset = gauss_function(normal_dist_energies[0], *popt)
			vertical_offset = 0
			for energy_val in normal_dist_energies:
				total_counts_from_gaussian.append(gauss_function(energy_val, *popt) - vertical_offset)
			angles_energies.append([i,peak_energy,sum(total_counts_from_gaussian)])
			#plt.text(peak_energy+50, peak_bins+100, r"$f\left(x,\mu,\sigma\right)=\frac{1}{\sigma\sqrt{2\pi}}e^{-\frac{\left(x-2\right)^{2}}{2\sigma^{2}}}$", fontsize=12)
			plt.text(170, (y_top)/2, r"$f\left(x,\mu,\sigma\right)=\frac{1}{"+str(sigma_neat)+r"\sqrt{2\pi}}e^{-\frac{\left(x-"+str(mean_neat)+r"\right)^{2}}{2\times"+str(sigma_neat)+r"^{2}}}$", fontsize=12)
			plt.text(670,(y_top)/2, r"$\int_{"+str(int(normal_dist_energies[0]))+r"}^{"+str(int(normal_dist_energies[len(normal_dist_energies)-1]))+r"}f\approx" + str(int(sum(total_counts_from_gaussian))) +r"$", fontsize="10")
			plt.plot(true_energy_vals, true_bin_vals, alpha=0.5, color="blue", linewidth="0.5", linestyle="-")
			plt.plot(normal_dist_energies,total_counts_from_gaussian, color="blue", linewidth="1")
			#plt.plot([peak_energy],[peak_bins],".", ms=2, color="blue", alpha=0.7)

		plt.savefig("all_subtracted_plots.png")
		plt.savefig("all_subtracted_plots.svg")
		plt.close()
		return angles_energies

	def get_angles_energies_plot(angles_energies_list, base_output_filename="angles_energies"):
		angles_list = []
		energies_list = []
		counts_list = []
		for coordinate in angles_energies_list:
			angles_list.append(coordinate[0])
			energies_list.append(coordinate[1])
			counts_list.append(coordinate[2])
		fit = np.polyfit(angles_list,energies_list,1)
		fit_fn = np.poly1d(fit)
		slope, intercept, r_value, p_value, std_err = linregress(angles_list,energies_list)
		plt.plot(angles_list,energies_list, "ro", ms=3)
		plt.plot(angles_list, fit_fn(angles_list), "--k")
		plt.ylabel("Energy (KeV)")
		plt.xlabel("Angle")
		plt.text(max(angles_list)/2,max(fit_fn(angles_list))/2, r"$E\prime="+("%0.2f"%slope)+r"\theta+"+("%0.2f"%intercept)+r"; R^{2} = "+("%0.2f"%(r_value**2)) + r"$", fontsize=12)
		plt.savefig(base_output_filename + ".png")
		plt.savefig(base_output_filename + ".svg")
		plt.close()
		#Now plot 1/E vs 1-cos(theta)
		cos_theta_angles = []
		for angle in angles_list:
			cos_theta_angles.append(1-math.cos(math.radians(angle)))
		inverse_energies_list_SI = []
		for energy in energies_list:
			new_energy = (1/float(energy))*1000
			inverse_energies_list_SI.append(new_energy)
		fit = np.polyfit(cos_theta_angles,inverse_energies_list_SI,1)
		fit_fn = np.poly1d(fit)
		slope, intercept, r_value, p_value, std_err = linregress(cos_theta_angles,inverse_energies_list_SI)
		#Electron eV is 1.782662 - find the % difference
		standard_ev = 1.782662
		percent_diff = abs((float(standard_ev)-float(slope))/float(standard_ev))*100
		ax = plt.gca()
		plt.text(0.5,0.25, r"$1/E^{\prime}="+("%0.2f"%slope)+r"\left(1-\cos\theta\right)+"+("%0.2f"%intercept)+r", R^{2} = "+("%0.2f"%(r_value**2)) + r"$", fontsize=12, transform=ax.transAxes)
		plt.text(0.5,0.15, "% Diff from expected electron mass: "+"{0:.2f}%".format(percent_diff), transform=ax.transAxes)
		plt.plot(cos_theta_angles, inverse_energies_list_SI, "ro", ms=3)
		plt.plot(cos_theta_angles, fit_fn(cos_theta_angles), "--k")
		plt.xlabel(r"$1-\cos\theta$", fontsize=12)
		plt.ylabel(r"$1/E^{\prime}$", fontsize=12)
		plt.title('Energy of Compton scattered gamma rays')
		plt.savefig("cos_theta_angles.svg")
		plt.savefig("cos_theta_angles.png")
		plt.close()
		#Now plot dsigma/dOmega vs. scattering angle theta
		differential_cross_section_values = []
		detector_efficiency_index = 0
		for count in counts_list:
			#Why do I need a factor of 2 to make my data look right?
			differential_cross_section_values.append(2*(count/detector_efficiencies[detector_efficiency_index]/720)/(1.0734)*(10**(-27))/100)
			detector_efficiency_index+=1
		angles_0_to_180 = np.arange(0,180,1)
		klein_nishina_values = []
		thompson_scattering_values = []
		c=299792458
		for theta in angles_0_to_180:
			klein_nishina_values.append(klein_nishina_function(0.6617, theta)*100)
			thompson_scattering_values.append(klein_nishina_function(0.6617/(c**2), theta)*100) #Klein Nishina fct converges to Thompson if E_gamma == m_e/c**2 http://en.wikipedia.org/wiki/Klein%E2%80%93Nishina_formula
		plt.ylabel(r"$\frac{d\sigma}{d\Omega}\mathrm{\left(10^{-28}cm^{2}/sr\right)}$", fontsize=12)
		plt.xlabel(r"$\mathrm{Scattering\,Angle\,\theta}$",fontsize=12)
		#plt.ylim(0,max[thompson_scattering_values])
		plt.plot(angles_0_to_180,klein_nishina_values)
		plt.plot(angles_0_to_180,thompson_scattering_values)
		plt.plot(angles_list,differential_cross_section_values, "bo", ms=3)
		#plt.scatter(angles_list,differential_cross_section_values, s=2, c="r")
		plt.savefig("differential_cross_section_plot.svg")
		plt.savefig("differential_cross_section_plot.png")
		plt.close()
	detector_efficiencies = [0.865,0.890, 0.930,0.945, 0.960, 0.975, 0.990, 0.9945,0.999]
	angles_energies = get_subtracted_plots()
	get_angles_energies_plot(angles_energies)

	get_all_plots()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Reads Quantum MCA data from Compton Scattering and generates relevant plots.')
	#parser.add_argument("-i", "--input")
	#parser.add_argument("-o", "--output")
	parser.add_argument("-p", "--peak-width-sample-size", default=50, help="Determines the sample size used for determining the width of the gaussian peak for each set of data. Default 50; the larger this value is, the wider the area of integration for gaussians will be")
	parser.add_argument("-a", "--averaging-threshold", default=28, help="For data that is negligible (i.e. around 0), when taking the average value over a sample of 100-200 data points, the data throughout that sampled interval will be considered as zero energy data. Default 28")
	args_dict = vars(parser.parse_args())
	main(args_dict)