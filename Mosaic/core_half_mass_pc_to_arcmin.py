import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate the projected core and half-mass radius on the sky in arcminutes.')
    parser.add_argument('-D', '--D_pc', type=float, required=True, help='Distance to the cluster in parsecs')
    parser.add_argument('-C', '--core_radius_pc', type=float, required=True, help='Core radius in parsecs')
    parser.add_argument('-H', '--half_mass_radius_pc', type=float, required=True, help='Half-mass radius in parsecs')
    return parser.parse_args()

def calculate_theta(radius_pc, D_pc):
    rad_to_arcmin = 3437.75
    pc_to_m = 3.09e16

    # Convert to meters
    radius_m = radius_pc * pc_to_m
    D_m = D_pc * pc_to_m

    theta_rad = radius_m / D_m
    theta_arcmin = theta_rad * rad_to_arcmin
    return theta_arcmin

def main():
    args = parse_arguments()

    # Calculate thetas
    theta_core_arcmin = calculate_theta(args.core_radius_pc, args.D_pc)
    theta_half_mass_arcmin = calculate_theta(args.half_mass_radius_pc, args.D_pc)

    print(f"The core radius is approximately {theta_core_arcmin:.6f} arcminutes.")
    print(f"The half-mass radius is approximately {theta_half_mass_arcmin:.6f} arcminutes.")

if __name__ == "__main__":
    main()
