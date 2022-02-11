# -*- coding: utf-8 -*-
import numpy as np
import pyik.adst as adst
from pyik.adst import RecEventProvider


input_file = '/home/bhanson/ADST_for_Ben.root'
output_folder = '/home/bhanson/Example_Events/'
for ev in RecEventProvider(input_file):
    event_id = str(ev.GetSDEvent().GetEventId())
    try:
        """ Read-out e.g. shower parameters """
        log10_rec_energy = np.log10(ev.GetSDEvent().GetSdRecShower().GetEnergy())
        log10_rec_energy_err = np.log10(ev.GetSDEvent().GetSdRecShower().GetEnergyError())
        log10_mc_energy = np.log10(ev.GetGenShower().GetEnergy())
        rec_shower_axis = (ev.GetSDEvent().GetSdRecShower().GetZenith(),
                       ev.GetSDEvent().GetSdRecShower().GetAzimuth())
        mc_shower_axis = (ev.GetGenShower().GetZenith(), 
                        ev.GetGenShower().GetAzimuth())
        mc_zenith, mc_azimuth = mc_shower_axis[0], mc_shower_axis[1]
        rec_zenith, rec_azimuth = rec_shower_axis[0], rec_shower_axis[1]
        mc_shower_core = np.asarray(ev.GetGenShower().GetCoreSiteCS())
        rec_shower_core = np.asarray(ev.GetSDEvent().GetSdRecShower().GetCoreSiteCS())
        trigger_status_vector = ev.GetSDEvent().GetStationVector()
        s1000 = ev.GetSDEvent().GetSdRecShower().GetS1000()
        s1000_err = ev.GetSDEvent().GetSdRecShower().GetS1000Error()


        dense_vector = []
        candidate_vector = []
        accident_vector = []
        rejected_vector = []
        silent_vector = []
        unknown_vector = []
        random_vector = []
        full_station_IDs = []
        tot_vector = []
        t1_vector = []
        t2_vector = []
        for station in trigger_status_vector:
            full_station_IDs.append(station.GetId())
            if station.IsDense():
                dense_vector.append(station.GetId())
            if station.IsCandidate():
                candidate_vector.append(station.GetId())
            if station.IsAccidental():
                accident_vector.append(station.GetId())
            if station.IsRejected():
                rejected_vector.append(station.GetId())
            if station.IsSilent():
                silent_vector.append(station.GetId())
            if station.IsUnknown():
                unknown_vector.append(station.GetId())
            if station.IsRandom():
                random_vector.append(station.GetId())
            if station.IsToT():
                tot_vector.append(station.GetId())
            if station.IsT1Threshold():
                t1_vector.append(station.GetId())
            if station.IsT2Threshold():
                t2_vector.append(station.GetId())           

                
        np.savez(output_folder + 'shower_parameters_'+ event_id + '.npz', event_id=int(event_id),
                 MC_Log_Energy=log10_mc_energy, Reconstructed_Log_Energy=log10_rec_energy,
                 Reconstructed_Log_Energy_Err = log10_rec_energy_err,
                 MC_Zenith = mc_zenith, Reconstructed_Zenith=rec_zenith, 
                 MC_Azimuth = mc_azimuth, Reconstructed_Azimuth=rec_azimuth, 
                 MC_Shower_Axis = mc_shower_axis, Reconstructed_Shower_Axis=rec_shower_axis,
                 MC_Shower_Core = mc_shower_core, Reconstructed_Shower_Core=rec_shower_core,
                 Trigger_Status_Count = len(np.asarray(trigger_status_vector)),
                 Full_Station_ID_Vector = np.asarray(full_station_IDs),
                 Candidate_Vector = np.asarray(candidate_vector), 
                 Dense_Vector = np.asarray(dense_vector), 
                 Accident_Vector = np.asarray(accident_vector), 
                 Rejected_Vector = np.asarray(rejected_vector), 
                 Random_Vector = np.asarray(random_vector), 
                 Silent_Vector = np.asarray(silent_vector), 
                 Unknown_Vector = np.asarray(unknown_vector),
                 ToThreshold_Vector = np.asarray(tot_vector), T1Threshold_Vector = np.asarray(t1_vector), T2Threshold_Vector = np.asarray(t2_vector),
                 S1000 = s1000, S1000_Err = s1000_err)
    except Exception as e:
        print('Something went wrong with event ' + event_id, str(e))
