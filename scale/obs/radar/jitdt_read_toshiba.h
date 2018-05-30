int jitdt_read_toshiba(int n_type, char *jitdt_place, pawr_header hd[n_type],
                       float az[n_type][ELDIM][AZDIM], float el[n_type][ELDIM][AZDIM],
                       float rtdat[n_type][ELDIM][AZDIM][RDIM]);
