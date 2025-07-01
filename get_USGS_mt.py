import pandas as pd

from obspy import UTCDateTime, Catalog
from obspy.clients.fdsn import Client

# =========================================================================

# Written by Amanda McPherson, June 2025
# Script to download all USGS moment tensor solutions within a lon/lat box
# of Alaska, sort the events, and output a .tsv file

# This can be used for any lon/lat box, time range, magnitude range, depth range, etc.

# Output columns:
# otime lon lat dep Mw(r or w) mag_type mrr mtt mff mrt mrf mtf contributor

# If there is not a Mwr or Mww magnitude, that column is filled with '--'
# =========================================================================

# File name to save results under
save_path = './USGS_Alaska_Catalog.tsv' # tab separator is hard-coded in

# Time range to search for events over

start = UTCDateTime('1999-12-31 23:59:59') # USGS Catalog starts around 2000
end = UTCDateTime() # Today's date

# Longitude and latitude box to search over
ax = [188, 228, 52, 72] # [lon_min, lon_max, lat_min, lat_max]


# Optional things to change

# Depth range to search over
dep_min = 0 # Can be None or a float in km
dep_max = 700 # Can be None or a float in km

# Magnitude range to search over
mag_min = 0.0 # Can be None
mag_max = 10.0 # Can be None


# Start searching for events

c = Client('USGS',timeout=600) # timeout is in seconds, default is 120

# Gets all events in the catalog within the specified ranges
# If doing the whole catalog, this will take a long time
# Once the number of events exceeds 20000, you will get an error and will have the batch download events
kwargs = {
    'starttime': start,
    'endtime': end,
    'minlatitude': ax[2],
    'maxlatitude': ax[3],
    'minlongitude': ax[0],
    'maxlongitude': ax[1],
    'mindepth' : dep_min,
    'maxdepth': dep_max,
    'minmagnitude': mag_min,
    'maxmagnitude': mag_max,
    'includeallorigins': True,
    'producttype': 'moment-tensor'
    }

print('Getting events...')
catalog = c.get_events(**kwargs)

print('Number of Events in catalog: ', len(catalog))

# Separate out the catalog information to put into a pandas DataFrame

lons, lats, deps, mags, mag_types = [], [], [], [], []
mrr, mtt, mff, mrt, mrf, mtf = [], [], [], [], [], []
otime = []
contrib = []

# USGS often has multiple solutions
# Hierarchy:
#    - Always tagged preferred solution
#    - USGS NEIC
#    - AEC
#    - Anything else

print('Organizing catalog information...')
for event in catalog:
    
    if event.preferred_focal_mechanism().moment_tensor:
    
        meth_id = event.preferred_focal_mechanism().moment_tensor.method_id
        
        try:
            mag_id = event.preferred_focal_mechanism().moment_tensor.moment_magnitude_id
        except:
            pass
        
        origin = event.preferred_origin().time
        otime.append(origin.datetime)
        lons.append(event.preferred_origin().longitude)
        lats.append(event.preferred_origin().latitude)
        deps.append(event.preferred_origin().depth / 1e3) # Depth in km

        # Specfically grab the magnitude associated with the preferred moment tensor if it is there
        if mag_id:
            all_mags = event.magnitudes
            for mag in all_mags:
                if mag.resource_id == mag_id:
                    mw = mag

            mag_type = mw.magnitude_type
            magnitude = mw.mag

        else:
            mag_type = event.preferred_magnitude().magnitude_type
            magnitude = event.preferred_magnitude().mag

            
        if mag_type.lower() in ['mwr', 'mww', 'mwb', 'mwc']:
            mags.append(magnitude)
            mag_types.append(mag_type)
        elif meth_id is not None:
            mags.append('--')
            mag_types.append(meth_id)
        else:
            mags.append('--')
            mag_types.append('--')

        fm = dict(event.preferred_focal_mechanism().moment_tensor.tensor)
        mrr.append(fm['m_rr'])
        mtt.append(fm['m_tt'])
        mff.append(fm['m_pp'])
        mrt.append(fm['m_rt'])
        mrf.append(fm['m_rp'])
        mtf.append(fm['m_tp'])

        contrib.append(event.preferred_focal_mechanism().creation_info.agency_id)
        
    else:
        fms = event.focal_mechanisms
        
        # Build a list of agency_ids
        ag_ids = []
        for ii in range(len(fms)):
            ag_ids.append(fms[ii].creation_info.agency_id)

        # Prefer NEIC, then AEC, then any solution, skip if can't find any solution
        try:
            neic = ag_ids.index('us')
            if fms[neic].moment_tensor.tensor:
                jj = neic
        except:
            try:
                aec = ag_ids.index('ak')
                if fms[aex].moment_tensor.tensor:
                    jj = aec
            except:
                for ii in range(len(fms)):
                    fm = fms[ii]
                    try:
                        if fm.moment_tensor.tensor:
                            jj = ii
                            break
                    except:
                        continue
    
            
        # Grab all the same information as above
        origin = event.preferred_origin().time
        otime.append(origin.datetime)
        lons.append(event.preferred_origin().longitude)
        lats.append(event.preferred_origin().latitude)
        deps.append(event.preferred_origin().depth / 1e3) # Depth in km
        
            
        # Specfically grab the magnitude associated with the moment tensor solution
        meth_id = fms[jj].moment_tensor.method_id
        try:
            mag_id = event.moment_tensor.moment_magnitude_id
        except:
            mag_id = False
            pass
        
        if mag_id:
            all_mags = event.magnitudes
            for mag in all_mags:
                if mag.resource_id == mag_id:
                    mw = mag

            mag_type = mw.magnitude_type
            magnitude = mw.mag

        else:
            mag_type = event.preferred_magnitude().magnitude_type
            magnitude = event.preferred_magnitude().mag
            
        if mag_type.lower() in ['mwr', 'mww', 'mwb', 'mwc']:
            mags.append(magnitude)
            mag_types.append(mag_type)
        elif meth_id is not None:
            mags.append('--')
            mag_types.append(meth_id)
        else:
            mags.append('--')
            mag_types.append('--')


        fm = dict(fms[jj].moment_tensor.tensor)
        mrr.append(fm['m_rr'])
        mtt.append(fm['m_tt'])
        mff.append(fm['m_pp'])
        mrt.append(fm['m_rt'])
        mrf.append(fm['m_rp'])
        mtf.append(fm['m_tp'])
    
        contrib.append(fms[jj].creation_info.agency_id)
            
    
# Organize everything into the desired organization, put into a pandas DataFrame, and save it as a .tsv
print('Final number of events: ' + str(len(otime)))
df = pd.DataFrame(list(zip(otime,lons,lats,deps,mags,mag_types,mrr,mtt,mff,mrt,mrf,mtf,contrib)))
print('Saving to .tsv...')
df.to_csv(save_path, sep='\t', index=False, header=False)
print('Finished.')
