import sys

def rotate_dump(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines = f_in.readlines()
        
        # State tracking
        in_atoms_section = False
        box_lines = []
        
        for line in lines:
            # 1. Handle Box Dimensions
            if "lo xhi" in line or "lo yhi" in line or "lo zhi" in line:
                box_lines.append(line.split())
                continue
            
            # Write out metadata and header
            if not in_atoms_section:
                f_out.write(line)
                
            # Detect Start of Atoms
            if "Atoms" in line:
                in_atoms_section = True
                
                # Process box dimensions for the swap
                # Index 0: lo, 1: hi, 2: label1, 3: label2
                x_box, y_box, z_box = box_lines[0], box_lines[1], box_lines[2]
                
                # Swap Y and Z bounds
                f_out.write(f"{x_box[0]} {x_box[1]} xlo xhi\n")
                f_out.write(f"{z_box[0]} {z_box[1]} ylo yhi\n")
                f_out.write(f"{y_box[0]} {y_box[1]} zlo zhi\n")
                continue
            
            # 2. Handle Atomic Coordinates
            if in_atoms_section and len(line.split()) >= 5:
                parts = line.split()
                # Atom structure: ID, Type, x, y, z
                # We want: ID, Type, x, z, y
                new_line = f"{parts[0]} {parts[1]} {parts[2]} {parts[4]} {parts[3]}\n"
                f_out.write(new_line)

if __name__ == "__main__":
    rotate_dump('./GaN_ortho.lmp', './GaN_ortho_rotated.lmp')
    print("Transformation complete: y and z columns swapped.")
