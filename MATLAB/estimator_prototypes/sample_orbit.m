[truth,measures] = sample_model(15*24*60*60*10,@orbit_init,@orbit_updater,@gps_measurer,0.1,0);