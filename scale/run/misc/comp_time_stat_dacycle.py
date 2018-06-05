import os


member = 100
group = member + 2
nnodes = 720


def parse_log_file(file):

    logtime = {}
    with open(file) as f:
        for line in f:
            if " *** ID=" in line:
                id = int(line[8:11])
                lab = line[14:47].strip()
                ave = float(line[54:64])
                max = float(line[73:83])
                min = float(line[99:109])
                logtime[lab] = {'id': id, 'ave': ave, 'max': max, 'min': min}

    return logtime


def scale_log_stat(logtime):

    stat = {}

    stat['init'] = logtime['INIT_Initialize']['max']
    stat['init_I'] = logtime['INIT_FILE_I_NetCDF']['max']
    stat['init_O'] = logtime['INIT_FILE_O_NetCDF']['max']
    stat['init_comm'] = logtime['INIT_COMM_vars']['max'] + logtime['INIT_COMM_wait']['max'] + logtime['INIT_COMM_Bcast']['max'] \
                      + logtime['INIT_COMM_init_pers']['max']
    stat['init_other'] = stat['init'] - stat['init_I'] - stat['init_O'] - stat['init_comm']

    stat['main'] = logtime['MAIN_Main_Loop']['max']
    stat['main_I'] = logtime['MAIN_FILE_I_NetCDF']['max']
    stat['main_O'] = logtime['MAIN_FILE_O_NetCDF']['max']
    stat['main_comm'] = logtime['MAIN_COMM_vars']['max'] + logtime['MAIN_COMM_wait']['max'] + logtime['MAIN_COMM_Allreduce']['max'] \
                      + logtime['MAIN_COMM_vars_pers']['max'] + logtime['MAIN_COMM_wait_pers']['max']
    stat['main_other'] = stat['main'] - stat['main_I'] - stat['main_O'] - stat['main_comm']

    return stat


def parse_NOUT_file(file):

    logtime = {}
    with open(file) as f:
        for line in f:
            if "##### TIMER" in line:
                id = line[14:64].strip()
                tt = float(line[78:92])
                if id in logtime:
                    logtime[id] += tt
                else:
                    logtime[id] = tt
            elif "#### TIMER" in line:
                id = line[17:67].strip()
                if (id.startswith('obsope_cal:read_ens_history')):
                    id = 'obsope_cal:read_ens_history'
                tt = float(line[81:95])
                if id in logtime:
                    logtime[id] += tt
                else:
                    logtime[id] = tt
            elif "### TIMER" in line:
                id = line[20:70].strip()
                if (id.startswith('write_restart_')):
                    continue
                tt = float(line[84:98])
                if id in logtime:
                    logtime[id] += tt
                else:
                    logtime[id] = tt

    return logtime


if __name__ == '__main__':

    import copy

#    import sys
#    no = sys.argv[1]

#    logpath = "20130713060000/log_{0:s}/scale".format(no)
    logpath = "20130713060000/log/scale"

    files = []
    for file in os.listdir(logpath):
        if file.endswith("_LOG.pe000000"):
            files += [file]

    files.sort()

    scale_stat = []
    m = 0
    mtot = 0
    g = 0
    for file in files:
        if m == 0:
            g += 1
            print "Group:", g
            time_total_max = 0.

        print "  File:", file

        logtime = parse_log_file(logpath + '/' + file)

        time_total = logtime['INIT_Initialize']['max'] + logtime['MAIN_Main_Loop']['max']
        if time_total > time_total_max:
            time_total_max = time_total
            logtime_max = logtime
            file_max = file

        m += 1
        mtot += 1

        if m == group or mtot == len(files):
            m = 0
            print 'Maximum total time:', time_total_max, '[Member', file_max[0:4], ']'
            stat_tmp = scale_log_stat(logtime_max)
            scale_stat += [stat_tmp]

            if g == 1:
                scale_stat_combine = copy.copy(stat_tmp)
            else:
                for key in stat_tmp.keys():
                    scale_stat_combine[key] += stat_tmp[key]


#    print ''
##    logpath = "20130713060000/log_{0:s}/scale".format(no)
#    logpath = "20130713060000/log/scale"

#    g = 0
#    gm = 0
#    for m in range(member+2):
##        filename = logpath + '/NOUT-' + str(g+1) + '.' + str(gm*nnodes)
#        filename = logpath + '/NOUT.' + str(gm*nnodes)
#        if gm == 0:
#            if os.path.isfile(filename):
#                print filename

#                logtime = parse_NOUT_file(filename)

#                if g == 0:
#                    scale_log_c = copy.copy(logtime)
#                else:
#                    for key in logtime.keys():
#                        scale_log_c[key] += logtime[key]

#        gm += 1
#        if gm >= group:
#            g += 1
#            gm = 0


#    print ''
#    logpath = "20130713060030/log/obsope"

#    g = 0
#    gm = 0
#    for m in range(member+2):
#        filename = logpath + '/NOUT.' + str(gm*nnodes)
#        if m == 0:
#            if os.path.isfile(filename):
#                print filename
#                obsope_log_m1 = parse_NOUT_file(filename)

#        gm += 1
#        if gm >= group:
#            g += 1
#            gm = 0


    print ''
#    logpath = "20130713060030/log_{0:s}/letkf".format(no)
    logpath = "20130713060000/log/dacycle"

    g = 0
    gm = 0
    for m in range(member+2):
        filename = logpath + '/NOUT.' + str(gm*nnodes)
        if m == 0:
            if os.path.isfile(filename):
                print filename
                letkf_log_m1 = parse_NOUT_file(filename)

        elif m == member:
            if os.path.isfile(filename):
                print filename
                letkf_log_mmean = parse_NOUT_file(filename)

        gm += 1
        if gm >= group:
            g += 1
            gm = 0


    print ''

#    items = ['PRE_SCRIPT', 'INITIALIZE', 'SCALE_RM', 'POST_SCRIPT']
#    for item in items:
#        print "{0:30s} {1:10.2f}".format(item, scale_log_c[item])

    items = ['init', 'init_I', 'init_O', 'init_comm', 'init_other', 'main', 'main_I', 'main_O', 'main_comm', 'main_other']
    for item in items:
        print "{0:30s} {1:10.2f}".format(item, scale_stat_combine[item])

    print ''

#    items = ['PRE_SCRIPT', 'INITIALIZE', 'READ_OBS']
#    for item in items:
#        print "{0:30s} {1:10.2f}".format(item, obsope_log_m1[item])

#    print "{0:30s} {1:10.2f}".format('read_history', obsope_log_m1['obsope_cal:read_ens_history_iter'])
#    print "{0:30s} {1:10.2f}".format('mean_calc', obsope_log_m1['OBS_OPERATOR'] - obsope_log_m1['obsope_cal:read_ens_history_iter'])

#    items = ['FINALIZE', 'POST_SCRIPT']
#    for item in items:
#        print "{0:30s} {1:10.2f}".format(item, obsope_log_m1[item])

#    print ''

    items = ['INITIALIZE', 'INIT_LETKF']
    for item in items:
        print "{0:30s} {1:10.2f}".format(item, letkf_log_m1[item])

    print "{0:30s} {1:10.2f}".format('SET_GRID(read_topo)', letkf_log_mmean['set_common_mpi_grid:read_topo:'])
    print "{0:30s} {1:10.2f}".format('SET_GRID(cal/scatter)', letkf_log_m1['SET_GRID'] - letkf_log_mmean['set_common_mpi_grid:read_topo:'])

    print ''

    print "{0:30s} {1:10.2f}".format('SCALE', letkf_log_m1['SCALE'])

    print "{0:30s} {1:10.2f}".format('READ_OBS(read)', letkf_log_m1['READ_OBS'] - letkf_log_m1['read_obs_all_mpi:mpi_bcast:'])
    print "{0:30s} {1:10.2f}".format('READ_OBS(bcast)', letkf_log_m1['read_obs_all_mpi:mpi_bcast:'])

    print "{0:30s} {1:10.2f}".format('OBS_OPERATOR(read_hist)', letkf_log_m1['obsope_cal:read_ens_history'])
    print "{0:30s} {1:10.2f}".format('OBS_OPERATOR(cal)', letkf_log_m1['OBS_OPERATOR'] - letkf_log_m1['obsope_cal:read_ens_history'])

    print "{0:30s} {1:10.2f}".format('PROCESS_OBS', letkf_log_m1['PROCESS_OBS'])

    print "{0:30s} {1:10.2f}".format('READ_GUES(read_rst)', letkf_log_m1['read_ens_mpi:read_restart:'])
    print "{0:30s} {1:10.2f}".format('READ_GUES(cal/scatter)', letkf_log_m1['READ_GUES'] - letkf_log_m1['read_ens_mpi:read_restart:'])

    print "{0:30s} {1:10.2f}".format('GUES_MEAN(cal/mon/gather)', letkf_log_m1['GUES_MEAN'] - letkf_log_mmean['write_ensmean:write_restart:'])
    print "{0:30s} {1:10.2f}".format('GUES_MEAN(write_rst)', letkf_log_m1['GUES_MEAN'])

    print "{0:30s} {1:10.2f}".format('DAS_LETKF', letkf_log_m1['DAS_LETKF'])

    print "{0:30s} {1:10.2f}".format('ANAL_MEAN(cal)', letkf_log_m1['ANAL_MEAN'])

    print "{0:30s} {1:10.2f}".format('WRITE_ANAL(write_rst)', letkf_log_m1['write_ens_mpi:write_restart:']) 
    print "{0:30s} {1:10.2f}".format('WRITE_ANAL(cal/mon/gather)', letkf_log_m1['WRITE_ANAL'] - letkf_log_m1['write_ens_mpi:write_restart:'])

#    print "{0:30s} {1:10.2f}".format('ensmspr_calc', letkf_log_mmean['GUES_MEAN'] + letkf_log_mmean['ANAL_MEAN'] -
#                                                     letkf_log_mmean['write_ensmspr_mpi:monit_obs'] - letkf_log_mmean['write_ensmspr_mpi:monit_obs_reduce_print'] -   
#                                                     letkf_log_mmean['write_ensmspr_mpi:write_mean'] - letkf_log_mmean['write_ensmspr_mpi:write_spread'])
#    print "{0:30s} {1:10.2f}".format('ensmspr_monit', letkf_log_mmean['write_ensmspr_mpi:monit_obs'] + letkf_log_mmean['write_ensmspr_mpi:monit_obs_reduce_print'])

    print ''

    items = ['FINALIZE',]
    for item in items:
        print "{0:30s} {1:10.2f}".format(item, letkf_log_m1[item])

