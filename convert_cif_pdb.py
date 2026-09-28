import os
from Bio.PDB import MMCIFParser, PDBIO

parser = MMCIFParser(QUIET=True)
io = PDBIO()

# loop through all files in current directory
for filename in os.listdir("."):
    if filename.endswith(".cif"):
        try:
            structure = parser.get_structure("structure", filename)
            
            output_filename = filename.replace(".cif", ".pdb")
            io.set_structure(structure)
            io.save(output_filename)
            
            print(f"✔ Converted: {filename} → {output_filename}")
        
        except Exception as e:
            print(f"✘ Failed: {filename} | Error: {e}")
