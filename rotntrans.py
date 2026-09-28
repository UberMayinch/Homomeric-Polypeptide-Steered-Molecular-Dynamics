from pymol import cmd
import os
import math
import re

# -----------------------------
# Config
# -----------------------------
PDB_FILES = [
    "decaglutamate.pdb",
    "decavaline.pdb",
    "decalysine.pdb",
    "decaglycine.pdb",
    "decaserine.pdb",
]

START_RESI = "1"
END_RESI = "10"
ATOM_NAME = "N"   # for N-to-N pulling. Change to "CA" if needed.

Z_LENGTH = 100.0  # arbitrary reference length along z-axis


def safe_object_name(filename):
    base = os.path.splitext(os.path.basename(filename))[0]
    return re.sub(r"[^A-Za-z0-9_]", "_", base)


def pick_atom(obj, resi, atom_name):
    """
    Prefer chain A if present; otherwise use residue number globally.
    """
    sel_chain_a = f"{obj} and chain A and resi {resi} and name {atom_name}"
    if cmd.count_atoms(sel_chain_a) == 1:
        return sel_chain_a

    sel_any_chain = f"{obj} and resi {resi} and name {atom_name}"
    count = cmd.count_atoms(sel_any_chain)

    if count != 1:
        raise RuntimeError(
            f"Expected exactly one atom for: {sel_any_chain}, found {count}"
        )

    return sel_any_chain


def check_z_alignment(start_sel, end_sel, label):
    a = cmd.get_atom_coords(start_sel)
    b = cmd.get_atom_coords(end_sel)

    dx = b[0] - a[0]
    dy = b[1] - a[1]
    dz = b[2] - a[2]

    length = math.sqrt(dx * dx + dy * dy + dz * dz)

    if length == 0:
        raise RuntimeError(f"{label}: start and end atoms overlap.")

    cos_theta = dz / length
    cos_theta = max(-1.0, min(1.0, cos_theta))
    angle_to_z = math.degrees(math.acos(cos_theta))

    print(f"\n{label}")
    print(f"  start coords: {a}")
    print(f"  end coords:   {b}")
    print(f"  vector: dx={dx:.6f}, dy={dy:.6f}, dz={dz:.6f}")
    print(f"  length: {length:.6f} Å")
    print(f"  angle to +z-axis: {angle_to_z:.6f} degrees")

    if abs(dx) < 1e-3 and abs(dy) < 1e-3 and dz > 0:
        print("  status: GOOD, aligned along +z")
    elif abs(dx) < 1e-3 and abs(dy) < 1e-3 and dz < 0:
        print("  status: aligned, but along -z")
    else:
        print("  status: NOT perfectly aligned")


def align_file_to_z(input_pdb):
    if not os.path.exists(input_pdb):
        print(f"Skipping missing file: {input_pdb}")
        return

    obj = safe_object_name(input_pdb)
    output_pdb = "aligned_" + os.path.basename(input_pdb)

    cmd.delete("all")
    cmd.load(input_pdb, obj)

    start_sel = pick_atom(obj, START_RESI, ATOM_NAME)
    end_sel = pick_atom(obj, END_RESI, ATOM_NAME)

    check_z_alignment(start_sel, end_sel, f"Before alignment: {input_pdb}")

    cmd.pseudoatom("zref", name="Z0", pos=[0.0, 0.0, 0.0])
    cmd.pseudoatom("zref", name="Z1", pos=[0.0, 0.0, Z_LENGTH])

    cmd.pair_fit(
        start_sel,
        "zref and name Z0",
        end_sel,
        "zref and name Z1",
    )

    cmd.delete("zref")

    start_sel = pick_atom(obj, START_RESI, ATOM_NAME)
    end_sel = pick_atom(obj, END_RESI, ATOM_NAME)

    check_z_alignment(start_sel, end_sel, f"After alignment: {output_pdb}")

    cmd.save(output_pdb, obj)
    print(f"  saved: {output_pdb}")


for pdb in PDB_FILES:
    align_file_to_z(pdb)

cmd.quit()
