import os
import shutil
import sys
import numpy as np


def generate_greens(args, data):
    delta = float(args.delta)
    number_of_cpus = float(args.cpus)
    number_of_samples = float(args.npts)
    reduction_velocity = float(args.rvel)

    # Load earth model
    # infile
    fin = open(os.path.join(args.info, "MODELS", args.model), 'r')
    lines = fin.readlines()
    fin.close()

    # read source earth model
    depths = []
    max_depth_model = 0
    for line in lines:
        tmp = line.rsplit(' ')
        max_depth_model += float(tmp[0])
        depths.append(float(tmp[0]))

    # max_depth_model -= float(tmp[0])

    # get depth value from args
    epi = args.epi.rsplit(' ')
    depth = float(epi[2])

    # if required depth too deep for model
    if depth > max_depth_model:
        print("Source too deep for your earth model. Exit!")
        sys.exit()

    # get layers from model
    number_of_layers = find_source_layer(depths, depth)

    # rewrite layers model with new depth
    number_of_layers, out_model = add_source_boundary(depths, lines, depth, list(number_of_layers))

    # get number of stations to generate Greens
    distances = get_the_distances(data)
    number_of_stations = distances.size
    # print(number_of_stations, distances)

    # Write first 13 lines of fkrprog Input file
    new_file_model = create_new_model(number_of_layers, out_model, depth, number_of_cpus, number_of_samples, delta,
                                      number_of_stations)

    # add station lines for each station to compute the Greens
    for i in range(number_of_stations):
        new_file_model.append("%8.4f%9.1f%10.1f" % (distances[i], 0, reduction_velocity))

    # Write new model file for fkrprog
    # out filename
    out_model_file = os.path.join(args.output, "output", 'MODEL1')
    fou = open(out_model_file, 'w')
    for i in range(len(new_file_model)):
        fou.write("%s\n" % (new_file_model[i]))
    fou.close()

    # remove old Greens file if exists and call fkprog to generate new greens
    # remove if exists
    if os.path.exists(os.path.join(args.output, "output", 'GREEN.1')):
        os.remove(os.path.join(args.output, "output", 'GREEN.1'))
    if os.path.exists('GREEN.1'):
        os.remove('GREEN.1')

    # commands to execute 1:cd outdoor 2: fkrprog
    # prog = os.path.join(args.info, "Programs", "FKRPROG.exe")
    cmd = "FKRPROG" + " < " + '"' + out_model_file + '"' + " > green.log"
    print(cmd)
    # subprocess.call([cmd])
    os.system(cmd)

    # move greens and log file
    if os.path.exists('green.log'):
        shutil.move('green.log', os.path.join(args.output, "output", 'green.log'))
    if os.path.exists('GREEN.1'):
        shutil.move('GREEN.1', os.path.join(args.output, "output", 'GREEN.1'))

    green = os.path.join(args.output, "output", 'GREEN.1')

    return green


def get_the_distances(data):
    distances = []

    for tr in data:
        distances.append(tr.stats['dist'])

    return np.array(list(set(distances)))


def find_source_layer(depths, depth):
    tmp = 0
    source_layer = 0
    while depth > tmp:
        tmp += depths[source_layer]
        source_layer += 1
    layers_below_source = len(depths) - source_layer
    return source_layer, layers_below_source


def add_source_boundary(depths, lines, depth, source_layer):
    # depths: array of depth from earth model (first column)
    # l: lines of models
    # depth: depth of source
    # source_layer: array of number of layers [0]: above the source and [1]: below the source

    for i in range(len(lines)):
        lines[i] = lines[i].rstrip()

    tmp = 0
    for i in range(source_layer[0]):
        tmp += float(depths[i])

    a = tmp - depth
    print(a)
    if abs(a) < 0.001:
        out_model = lines
        source_layer[1] = source_layer[1] + 0

    else:
        # reconstruct model
        # if b!=0
        # line 0 - nrs[0] ok
        # line nrs[0] split: top below with respect to the source depth
        # add nrs[0]-nrs[1] layers
        # update nrs[1]

        new_l = []
        tmp = 0
        for i in range(source_layer[0]):
            tmp += float(depths[i])

        for i in range(source_layer[0] - 1):
            new_l.append(lines[i])

        # extract depth from layer nrs[0]-1
        # and split into two layers according to depth source
        foe = lines[source_layer[0] - 1].rsplit(' ')

        # new thickness of layer
        thickness_below_depth = tmp - depth
        thickness_upper_depth = float(foe[0]) - thickness_below_depth

        layer_top_tmp = lines[source_layer[0] - 1]

        line_of_top = layer_top_tmp.rsplit(' ')
        line_of_bellow = layer_top_tmp.rsplit(' ')

        line_of_top[0] = "%.4E" % thickness_upper_depth
        line_of_bellow[0] = "%.4E" % thickness_below_depth

        line_of_top = ' '.join(line_of_top)
        line_of_bellow = ' '.join(line_of_bellow)

        new_l.append(line_of_top)
        new_l.append(line_of_bellow)

        for i in range(source_layer[1]):
            new_l.append(lines[source_layer[0] + i])

        out_model = new_l
        source_layer[1] = source_layer[1] + 1

    return source_layer, out_model


def create_new_model(number_of_layers, lines, depth, number_of_cpus, number_of_samples, delta, number_of_stations):
    # Station parameters
    ph_vel = [10000, 30]
    # ra_vel=[1.2,1.0]    # less than ralay veleocity of the model
    ra_vel = [2.6, 2.5]  # less than ralay veleocity of the model
    x = 6
    beg_freq = 1
    nr_l = number_of_layers[0] + number_of_layers[1]
    nr_fr = int(number_of_samples / 2)
    out = ['.F.', '    0   64', "GREEN.%1d" % number_of_cpus, "%7.1f%10.2f%8d%5d%5d%9.3f%6d%5d" % (
        x, depth, beg_freq, nr_fr, number_of_samples, delta, nr_l, number_of_cpus),
           "%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d" % (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)]

    # start writing model list
    # line 4
    # line 5
    # layers
    for i in range(len(lines)):
        out.append(" %s" % (lines[i]))
    out.append("%5d" % (number_of_layers[0] + 1))  # line M+1
    # Here to correct previous nrs[1] which was made for
    # "Line 11 gives the layer number below the source"
    # in the manual, but ...
    # "Line 11 gives the layer number above the source + 1 layer"
    out.append("%15.7E%14.6E%10d" % (400, 1.5, 0))  # line M+2
    # line M+3
    out.append("%5d%9.1f%9.1f%9.1f%10.1f" % (number_of_stations, ph_vel[0], ph_vel[1], ra_vel[0], ra_vel[1]))

    return out
