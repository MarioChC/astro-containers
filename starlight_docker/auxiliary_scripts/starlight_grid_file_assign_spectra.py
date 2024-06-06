import os
import argparse

def get_kinematic_info(input_file, kinematic_file):
    kinematic_info = {}
    with open(kinematic_file, 'r') as f_kinematic:
        kinematic_data = f_kinematic.readlines()
        for line in kinematic_data:
            line_data = line.split()
            if line_data:
                spectrum_name = line_data[0]
                if spectrum_name in input_file:
                    kinematic_info[spectrum_name] = (float(line_data[1]), float(line_data[3]))  # Store the velocity and sigma values
    return kinematic_info


def generate_starlight_input_file(output_file, spectrum_names, kinematic_file):

    with open(output_file, 'w') as f_output:
        f_output.write("2702                                             [Number of fits to run]\n")
        f_output.write("/home/starlight/STARLIGHTv04/                    [base_dir]\n")
        f_output.write("/home/starlight/shared_directory/config_files_starlight/spectrum/                        [obs_dir]\n")
        f_output.write("/home/starlight/shared_directory/config_files_starlight/mask/                            [mask_dir]\n")
        f_output.write("/home/starlight/shared_directory/output/                             [out_dir]\n")
        f_output.write("-2007200                                         [your phone number]\n")
        f_output.write("4730.0                                           [llow_SN]   lower-lambda of S/N window\n")
        f_output.write("4780.0                                           [lupp_SN]   upper-lambda of S/N window\n")
        f_output.write("2700.0                                           [Olsyn_ini] lower-lambda for fit\n")
        f_output.write("9000.0                                           [Olsyn_fin] upper-lambda for fit\n")
        f_output.write("1.0                                              [Odlsyn]    delta-lambda for fit\n")
        f_output.write("1.0                                              [fscale_chi2] fudge-factor for chi2\n")
        f_output.write("FXK                                              [FIT/FXK] Fit or Fix kinematics\n")
        f_output.write("0                                                [IsErrSpecAvailable]  1/0 = Yes/No\n")
        f_output.write("0                                                [IsFlagSpecAvailable] 1/0 = Yes/No\n")
        # Añadir nuevas líneas al final del archivo
        for spectrum_name_with_ext in spectrum_names:
            spectrum_name = os.path.splitext(os.path.basename(spectrum_name_with_ext))[0]  # Obtener solo el nombre del archivo sin la extensión ni la ruta
            kinematic_info = get_kinematic_info(spectrum_name_with_ext, kinematic_file)
            spectrum_name, (vel, sigma) = kinematic_info.popitem()
            new_line = f"{spectrum_name}.txt StCv04.C11.arp220.config Base.BC03.N Masks.Em.Abs.Lines.Arp220.gm CAL {vel} {sigma} {spectrum_name}_output.txt\n"
            f_output.write(new_line)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Modify a file and add new lines at the end")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("spectrum_names", nargs="*", help="Names of spectrum files")
    parser.add_argument("kinematic_file", help="Path to the file that contains the information about the kinematics")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    generate_starlight_input_file(args.output_file, args.spectrum_names, args.kinematic_file)