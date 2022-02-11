import numpy as np
import scipy.stats
import os, os.path
import matplotlib.pyplot as plt
import math
from collections import Counter
from scipy.stats import binned_statistic_2d
import itertools as IT
from scipy.stats import norm
from scipy import optimize

class Shower_Data:
    def __init__(self):
        self.event_id = []
        self.trigger_status = []
        self.mc_log_energy = []
        self.reconstructed_log_energy = []
        self.reconstructed_log_energy_err = []
        self.mc_zenith = []
        self.reconstructed_zenith= []
        self.mc_azimuth = []
        self.reconstructed_azimuth= []
        self.mc_shower_axis = []
        self.reconstructed_shower_axis = []
        self.mc_shower_core = []
        self.reconstructed_shower_core = []
        self.full_station_id_vector = []
        self.candidate_vector = []
        self.dense_vector = []
        self.accident_vector = []
        self.rejected_vector = []
        self.silent_vector = []
        self.unknown_vector = []
        self.random_vector = []
        self.tot_vector = []
        self.t1_vector = []
        self.t2_vector = []
        self.s1000 = []
        self.s1000_err = []

    def printAll(self):
        print("Event ID:")
        print(self.event_id)
        print("Trigger Status:")
        print(self.trigger_status)
        print("MC Log Energy")
        print(self.mc_log_energy)
        print("Reconstructed Log Energy")
        print(self.reconstructed_log_energy)
        print("Reconstructed Log Energy Error")
        print(self.reconstructed_log_energy_err)
        print("MC Zenith")
        print(self.mc_zenith)
        print("Reconstructed Zenith")
        print(self.reconstructed_zenith)
        print("MC Azimuth")
        print(self.mc_azimuth)
        print("Reconstructed Azimuth")
        print(self.reconstructed_azimuth)
        print("MC Shower Axis")
        print(self.mc_shower_axis)
        print("Reconstructed Shower Axis")
        print(self.reconstructed_shower_axis)
        print("MC Shower Core")
        print(self.mc_shower_core)
        print("Reconstructed Shower Core")
        print(self.reconstructed_shower_core)
        print("Full Station ID Vector")
        print(self.full_station_id_vector)
        print("Candidate Vector")
        print(self.candidate_vector)
        print("Dense Vector")
        print(self.dense_vector)
        print("Accident Vector")
        print(self.accident_vector)
        print("Rejected Vector")
        print(self.rejected_vector)
        print("Silent Vector")
        print(self.silent_vector)
        print("Unknown Vector")
        print(self.unknown_vector)
        print("Random Vector")
        print(self.random_vector)
        print("S1000")
        print(self.s1000)
        print("S1000 Error")
        print(self.s1000_err)
        print("ToThreshold")
        print(self.tot_vector)
        print("T1Threshold")
        print(self.t1_vector)
        print("T2Threshold")
        print(self.t2_vector)

    def printSingle(self, event_number):
        event_number = str(event_number)
        index = 0
        i = 0
        for i in range(len(self.event_id)):
            if (event_number == str(self.event_id[i])):
                index = i
                break
        print("Event ID:")
        print(self.event_id[i])
        print("Trigger Status:")
        print(self.trigger_status[i])
        print("MC Log Energy")
        print(self.mc_log_energy[i])
        print("Reconstructed Log Energy")
        print(self.reconstructed_log_energy[i])
        print("Reconstructed Log Energy Err")
        print(self.reconstructed_log_energy_err)
        print("MC Zenith")
        print(self.mc_zenith[i])
        print("Reconstructed Zenith")
        print(self.reconstructed_zenith[i])
        print("MC Azimuth")
        print(self.mc_azimuth[i])
        print("Reconstructed Azimuth")
        print(self.reconstructed_azimuth[i])
        print("MC Shower Axis")
        print(self.mc_shower_axis[i])
        print("Reconstructed Shower Axis")
        print(self.reconstructed_shower_axis[i])
        print("MC Shower Core")
        print(self.mc_shower_core[i])
        print("Reconstructed Shower Core")
        print(self.reconstructed_shower_core[i])
        print("Full Station ID Vector")
        print(self.full_station_id_vector[i])
        print("Candidate Vector")
        print(self.candidate_vector[i])
        print("Dense Vector")
        print("Accident Vector")
        print(self.accident_vector[i])
        print("Rejected Vector")
        print(self.rejected_vector[i])
        print("Silent Vector")
        print(self.silent_vector[i])
        print("Unknown Vector")
        print(self.unknown_vector[i])
        print("Random Vector")
        print(self.random_vector[i])
        print("S1000")
        print(self.s1000[i])
        print("S1000 Error")
        print(self.s1000_err[i])
        print("ToThreshold")
        print(self.tot_vector[i])
        print("T1Threshold")
        print(self.t1_vector[i])
        print("T2Threshold")
        print(self.t2_vector[i])

def make_plots(sorted_Data, plot_list):
    for i in range(len(plot_list)):
        if (plot_list[i] == 1):
            candidate_count = []
            for i in range(len(sorted_Data.candidate_vector)):
                candidate_count.append(len(sorted_Data.candidate_vector[i]))

            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Number of Candidate Stations vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Log of Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Number of Candidate Stations', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            xbins = np.linspace(18.5, 19.0, 50)
            ybins = np.linspace(0, 12, 13)
            mean_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, candidate_count, statistic = 'mean', bins = xbins)
            std_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, candidate_count, statistic = 'std', bins = xbins)
            means = mean_results.statistic
            std = std_results.statistic
            bin_edges = mean_results.bin_edges
            bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
            plt.errorbar(x = bin_centers, y = means, yerr = std, 
                        fmt = 'ro',
                        ms = 8,
                        ecolor = 'black',
                        elinewidth = 1.5,			
                        capsize = 3, 
                        capthick = 1.5, 
                        marker = '_') 
            plt.hist2d(sorted_Data.mc_log_energy, candidate_count, bins = [xbins, ybins], cmap='GnBu')
            plt.colorbar()
            plt.savefig('candidate_vs_energy.jpg', bbox_inches=0, dpi=600)
            plt.show()

        elif (plot_list[i] == 2):
            # Energy Reconstructed vs. Energy Thrown
            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Energy Reconstructed vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Log of Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Log of Reconstructed Energy (eV)', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.hist2d(sorted_Data.mc_log_energy,sorted_Data.reconstructed_log_energy, 50, cmap='GnBu')
            plt.colorbar()
            corrected_error = list(np.zeros(len(sorted_Data.reconstructed_log_energy_err)))
            for i in range(len(corrected_error)):
                corrected_error[i] = abs(sorted_Data.reconstructed_log_energy[i] - math.log10((10**sorted_Data.reconstructed_log_energy[i])+(10**sorted_Data.reconstructed_log_energy_err[i])))
            mean_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, sorted_Data.reconstructed_log_energy, statistic = 'mean', bins = 50)
            error_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, corrected_error, statistic = 'mean', bins = 50)
            means = mean_results.statistic
            error = error_results.statistic
            bin_edges = mean_results.bin_edges
            bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
            plt.errorbar(x = bin_centers, y = means, yerr = error, 
                        fmt = 'ro',
                        ms = 8,
                        ecolor = 'black',
                        elinewidth = 1.5,			
                        capsize = 3, 
                        capthick = 1.5, 
                        marker = '_') 
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            x = np.linspace(min(sorted_Data.mc_log_energy), max(sorted_Data.mc_log_energy), 10)
            y = np.linspace(min(sorted_Data.mc_log_energy), max(sorted_Data.mc_log_energy), 10)
            plt.plot(x, y, '--')
            plt.savefig('reconstructed_vs_thrown.jpg', bbox_inches=0, dpi=600)
            plt.show()
        elif (plot_list[i] == 3):
            # ToThreshold vs. Energy
            # Getting Total Stations
            total_stations = []
            for i in range(len(sorted_Data.full_station_id_vector)):
                total_stations.append([x for x in sorted_Data.full_station_id_vector[i] if x not in sorted_Data.dense_vector[i] and x not in sorted_Data.accident_vector[i] and x not in sorted_Data.unknown_vector[i] and x not in sorted_Data.rejected_vector[i] and x not in sorted_Data.silent_vector[i] and x not in sorted_Data.random_vector[i]])
                sorted_Data.tot_vector[i] = [x for x in sorted_Data.tot_vector[i] if x not in sorted_Data.dense_vector[i] and x not in sorted_Data.accident_vector[i] and x not in sorted_Data.unknown_vector[i] and x not in sorted_Data.rejected_vector[i] and x not in sorted_Data.silent_vector[i] and x not in sorted_Data.random_vector[i]]
                sorted_Data.t1_vector[i] = [x for x in sorted_Data.t1_vector[i] if x not in sorted_Data.dense_vector[i] and x not in sorted_Data.accident_vector[i] and x not in sorted_Data.unknown_vector[i] and x not in sorted_Data.rejected_vector[i] and x not in sorted_Data.silent_vector[i] and x not in sorted_Data.random_vector[i]]
                sorted_Data.t2_vector[i] = [x for x in sorted_Data.t2_vector[i] if x not in sorted_Data.dense_vector[i] and x not in sorted_Data.accident_vector[i] and x not in sorted_Data.unknown_vector[i] and x not in sorted_Data.rejected_vector[i] and x not in sorted_Data.silent_vector[i] and x not in sorted_Data.random_vector[i]]


            tot_percentage = []
            for i in range(len(total_stations)):
                tot_percentage.append((len(sorted_Data.tot_vector[i])/len(total_stations[i])*100))
            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Percentage of ToT Triggered Stations vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Log of Throw Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Percent of Stations ToT Triggered', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            xbins = np.linspace(18.5, 19.0, 50)
            ybins = np.linspace(0, 100, 20)
            mean_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, tot_percentage, statistic = 'mean', bins = xbins)
            std_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, tot_percentage, statistic = 'std', bins = xbins)
            means = mean_results.statistic
            std = std_results.statistic
            bin_edges = mean_results.bin_edges
            bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
            plt.errorbar(x = bin_centers, y = means, yerr = std, 
                        fmt = 'ro',
                        ms = 8,
                        ecolor = 'black',
                        elinewidth = 1.5,			
                        capsize = 3, 
                        capthick = 1.5, 
                        marker = '_') 
            plt.hist2d(sorted_Data.mc_log_energy, tot_percentage, bins = [xbins, ybins], cmap='GnBu')
            plt.colorbar()
            plt.savefig('tot_trigger_vs_energy.jpg', bbox_inches=0, dpi=600)
            plt.show()


            t1_percentage = []
            for i in range(len(total_stations)):
                t1_percentage.append((len(sorted_Data.t1_vector[i])/len(total_stations[i]))*100)
            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Percentage of T1 Triggers vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Log of Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Percentage of Stations T1 Triggered', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            xbins = np.linspace(18.5, 19.0, 50)
            ybins = np.linspace(0, 100, 20)
            mean_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, t1_percentage, statistic = 'mean', bins = xbins)
            std_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, t1_percentage, statistic = 'std', bins = xbins)
            means = mean_results.statistic
            std = std_results.statistic
            bin_edges = mean_results.bin_edges
            bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
            plt.errorbar(x = bin_centers, y = means, yerr = std, 
                        fmt = 'ro',
                        ms = 8,
                        ecolor = 'black',
                        elinewidth = 1.5,			
                        capsize = 3, 
                        capthick = 1.5, 
                        marker = '_') 
            plt.hist2d(sorted_Data.mc_log_energy,t1_percentage, bins = [xbins, ybins], cmap='GnBu')
            plt.colorbar()
            plt.savefig('t1_trigger_vs_energy.jpg', bbox_inches=0, dpi=600)
            plt.show()


            t2_percentage = []
            for i in range(len(total_stations)):
                t2_percentage.append((len(sorted_Data.t2_vector[i])/len(total_stations[i]))*100)
            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Percentage of T2 Triggers vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Log of Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Percentage of Stations T2 Triggered', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            xbins = np.linspace(18.5, 19.0, 50)
            ybins = np.linspace(0, 100, 20)
            mean_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, t2_percentage, statistic = 'mean', bins = xbins)
            std_results = scipy.stats.binned_statistic(sorted_Data.mc_log_energy, t2_percentage, statistic = 'std', bins = xbins)
            means = mean_results.statistic
            std = std_results.statistic
            bin_edges = mean_results.bin_edges
            bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
            plt.errorbar(x = bin_centers, y = means, yerr = std, 
                        fmt = 'ro',
                        ms = 8,
                        ecolor = 'black',
                        elinewidth = 1.5,			
                        capsize = 3, 
                        capthick = 1.5, 
                        marker = '_') 
            plt.hist2d(sorted_Data.mc_log_energy,t2_percentage, bins = [xbins, ybins], cmap='GnBu')
            plt.colorbar()
            plt.savefig('t2_trigger_vs_energy.jpg', bbox_inches=0, dpi=600)
            plt.show()
        elif (plot_list[i] == 4):
            #Zenith vs. Energy vs. Candidate Counts
            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Candidate Counts vs. Thrown Energy vs. Zenith Angle', fontsize = 25)
            plt.xlabel('Log of Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Cos^2(Zenith Angle (rad))', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            
            candidate_count = []
            for i in range(len(sorted_Data.candidate_vector)):
                candidate_count.append(len(sorted_Data.candidate_vector[i]))

            x = sorted_Data.mc_log_energy
            y = sorted_Data.mc_zenith
            z = candidate_count

            xbins = np.linspace(18.5, 19.0, 25)
            ybins = np.linspace(0.292, 1.0, 25)

            ret = binned_statistic_2d(x, y, z, statistic=np.mean, bins= [xbins, ybins])
            array = ret.statistic
            array = array.T
            array = np.nan_to_num(array)
            plt.imshow(array, origin='lower', extent = [18.5, 19.0, 0.292, 1.0], cmap = 'GnBu')
            plt.colorbar(label = 'Number of Candidate Stations')
            axes = plt.gca()
            axes.set_aspect(0.706)
            plt.savefig('candidate_counts_vs_zenith_vs_energy.jpg', bbox_inches=0, dpi=600)
            plt.show()
        elif (plot_list[i] == 5):
            #Making the Number of Counts vs. Energy
            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Number of Events vs. Thrown Energy vs. Zenith Angle', fontsize = 25)
            plt.xlabel('Log of Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Cos^2(Zenith Angle (rad))', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)

            x = sorted_Data.mc_log_energy
            y = list((np.cos(np.array(sorted_Data.mc_zenith)))**2)

            xbins = np.linspace(18.5, 19.0, 25)
            ybins = np.linspace(0.292, 1.0, 25)

            plt.hist2d(x,y, bins = [xbins, ybins], cmap = "GnBu")
            axes = plt.gca()
            axes.set_aspect(0.706)
            plt.colorbar()
            plt.savefig('events_vs_zenith_vs_energy.jpg', bbox_inches=0, dpi=600)
            plt.show()
        elif (plot_list[i] == 6):
            h, xedges, yedges, _ = plt.hist2d(sorted_Data.mc_log_energy,sorted_Data.reconstructed_log_energy, 10)
            plt.clf()

            mc_energies = []
            for i in range(len(xedges)):
                mc_energies.append(10**(xedges[i]))

            rec_energies = []
            for i in range(len(yedges)):
                rec_energies.append(10**(yedges[i]))

            mc_energies.pop()
            rec_energies.pop()

            resolution = []
            resolution_sub = []
            for i in range(len(mc_energies)):
                for j in range(len(rec_energies)):
                    resolution_sub.append((rec_energies[j] - mc_energies[i])/mc_energies[i])
                resolution.append(resolution_sub)
                resolution_sub = []


            sorted_h = []
            for i in range(len(h)):
                sorted_h.append([x for _, x in sorted(zip(resolution[i], h[i]))])
                resolution[i] = sorted(resolution[i])

            mean = []
            stddev = []
            fwhm = []


            fig=plt.figure()
            fig.suptitle('Resolution Distributions vs. Bins')
            for i in range(len(sorted_h)):
                x = resolution[i]
                y = sorted_h[i]
                popt, _ = optimize.curve_fit(gaussian, resolution[i], sorted_h[i])
                mean.append(popt[1])
                stddev.append(abs(popt[2]))
                fwhm.append(2*abs(popt[2])*0.83255)
                
                gauss = gaussian(x, *popt)

                ax = plt.subplot(2,5,i+1)
                ax.scatter(x,y, c = "blue", s = 2)
                ax.plot(x, gauss, c = "red")
            plt.savefig('resolution_vs_bins.jpg', bbox_inches=0, dpi=600)
            plt.show()

            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('FWHM vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('FWHM (resolution)', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.scatter(mc_energies, fwhm)
            plt.savefig('fwhm_vs_thrown.jpg', bbox_inches=0, dpi=600)
            plt.show()

            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Stddev vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Stddev (resolution)', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.scatter(mc_energies, stddev)
            plt.savefig('stddev_vs_thrown.jpg', bbox_inches=0, dpi=600)
            plt.show()

            fig= plt.figure(figsize=(12,8)) # Making plot bigger
            plt.title('Offset vs. Thrown Energy', fontsize = 25)
            plt.xlabel('Log of Thrown Energy (eV)', fontsize = 20) # Labeling the x-axis
            plt.ylabel('Offset/Mean (resolution)', fontsize = 20) # Labeling the y-axis # Specifying x-plot range
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.scatter(mc_energies, mean)
            plt.savefig('mean_vs_thrown.jpg', bbox_inches=0, dpi=600)
            plt.show()



def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)

sorted_Data = Shower_Data()
input_folder = 'C:/Users/19136/OneDrive - Colorado School of Mines/PHGN 471/Auger@TA/Example_Events_2/'
path = 'C:/Users/19136/OneDrive - Colorado School of Mines/PHGN 471/Auger@TA/Example_Written_Files/'

# Initialize Events - Only do this once
count = 0
for filename in os.listdir(input_folder):
    data = np.load(input_folder + filename)
    lst = data.files
    for item in lst:
        if (item == 'event_id'):
            sorted_Data.event_id.append(int(data[item]))
        elif (item == "Trigger_Status_Count"):
            sorted_Data.trigger_status.append(int(data[item]))
        elif (item == "MC_Log_Energy"):
            sorted_Data.mc_log_energy.append(round(float(data[item]), 2))
        elif (item == "Reconstructed_Log_Energy"):
            sorted_Data.reconstructed_log_energy.append(round(float(data[item]),2))
        elif (item == "Reconstructed_Log_Energy_Err"):
            sorted_Data.reconstructed_log_energy_err.append(round(float(data[item]),2))
        elif (item == "MC_Zenith"):
            sorted_Data.mc_zenith.append(round(float(data[item]),2))
        elif (item == "Reconstructed_Zenith"):
            sorted_Data.reconstructed_zenith.append(round(float(data[item]),2))
        elif (item == "MC_Azimuth"):
            sorted_Data.mc_azimuth.append(round(float(data[item]),2))
        elif (item == "Reconstructed_Azimuth"):
            sorted_Data.reconstructed_azimuth.append(round(float(data[item]),2))
        elif (item == "MC_Shower_Axis"):
            sorted_Data.mc_shower_axis.append(list(np.around(data[item],2)))
        elif (item == "Reconstructed_Shower_Axis"):
            sorted_Data.reconstructed_shower_axis.append(list(np.around(data[item],2)))
        elif (item == "MC_Shower_Core"):
            sorted_Data.mc_shower_core.append(list(np.around(data[item],2)))
        elif (item == "Reconstructed_Shower_Core"):
            sorted_Data.reconstructed_shower_core.append(list(np.around(data[item],2)))
        elif(item == "Full_Station_ID_Vector"):
            sorted_Data.full_station_id_vector.append(list(data[item]))
        elif(item == "Candidate_Vector"):
            sorted_Data.candidate_vector.append(list(data[item]))
        elif(item == "Dense_Vector"):
            sorted_Data.dense_vector.append(list(data[item]))
        elif(item == "Accident_Vector"):
            sorted_Data.accident_vector.append(list(data[item]))
        elif(item == "Rejected_Vector"):
            sorted_Data.rejected_vector.append(list(data[item]))
        elif (item == "Silent_Vector"):
            sorted_Data.silent_vector.append(list(data[item]))
        elif (item == "Unknown_Vector"):
            sorted_Data.unknown_vector.append(list(data[item]))
        elif (item == "Random_Vector"):
            sorted_Data.random_vector.append(list(data[item]))
        elif(item == "S1000"):
            sorted_Data.s1000.append(int(data[item]))
        elif(item == "S1000_Err"):
            sorted_Data.s1000_err.append(int(data[item]))
        elif(item == "ToThreshold_Vector"):
            sorted_Data.tot_vector.append(list(data[item]))
        elif(item == "T1Threshold_Vector"):
            sorted_Data.t1_vector.append(list(data[item]))
        elif(item == "T2Threshold_Vector"):
            sorted_Data.t2_vector.append(list(data[item]))
        else:
            print("Item identifier not found: ", item)

'''Plot types:
        1. Number of Candidate Stations vs. Energy
        2. Energy Reconstructed vs. Energy Thrown
        3. Trigger types vs. Energy
        4. Number of Candidate Stations vs. Zenith Angle vs. Thrown Energy 
        5. Number of Events vs. Zenith Angle vs. Thrown Energy
        6. Std Dev/FWHM/Offset Resolution vs. Thrown Energy'''


plot_list = [4,5]
make_plots(sorted_Data, plot_list)


