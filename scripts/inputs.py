import argparse

#Define arguments for each required and optional input
def define_arguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ## Required inputs
    parser.add_argument("--input-mech-file",dest="InputMechFile",required=True,help="InputMechFile")

    parser.add_argument("--labels-file",dest="LabelsFile",required=True,help="LabelsFile")

    return parser

# Wrapper function
def generate_arguments():

    #Generate argument parser and define arguments
    parser = define_arguments()
    args = parser.parse_args()
    input_mech_file = args.InputMechFile
    labels_file = args.LabelsFile

    return input_mech_file,labels_file


