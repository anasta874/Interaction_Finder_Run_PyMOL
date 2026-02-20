import os
import math
from pymol import cmd

# =========================
# CONFIG 
# =========================

# Output
output_directory = r"C:/path/to/output"
output_filename = "interactions_list.txt"

# Object names in PyMOL (MUST be separate objects)
rec_file = "rec"
lig_file = "lig"

# Cutoffs
distance_cutoff = 4.0
hbond_distance_cutoff = 3.5
hbond_angle_cutoff = 120.0

ensure_hydrogens = False


# HELPERS

def clamp(x, lo=-1.0, hi=1.0):
    return max(lo, min(hi, x))

def calculate_angle(atom1, atom2, atom3):
    """Angle (deg) for atom1-atom2-atom3 using coordinates."""
    v1 = (
        atom1.coord[0] - atom2.coord[0],
        atom1.coord[1] - atom2.coord[1],
        atom1.coord[2] - atom2.coord[2],
    )
    v2 = (
        atom3.coord[0] - atom2.coord[0],
        atom3.coord[1] - atom2.coord[1],
        atom3.coord[2] - atom2.coord[2],
    )

    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    if mag1 == 0 or mag2 == 0:
        return 0.0

    cosang = clamp(dot / (mag1 * mag2))
    return math.degrees(math.acos(cosang))

def distance_coords(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def get_element(atom):
    """Prefer atom.symbol; fallback to atom.name heuristics."""
    sym = getattr(atom, "symbol", None)
    if sym:
        return sym.strip().upper()

    name = (atom.name or "").strip().upper()
    # crude fallback for two-letter elements commonly seen in PDB
    if len(name) >= 2 and name[:2] in {"CL", "BR", "NA", "CA", "MG", "ZN", "FE", "MN", "CU", "NI", "CO"}:
        return name[:2]
    return name[:1] if name else ""

def ensure_object_exists(obj_name):
    objs = cmd.get_object_list()
    if obj_name not in objs:
        raise ValueError(f"Object '{obj_name}' not found in PyMOL. Loaded objects: {objs}")

def ensure_h_atoms(obj_name):
    """Add hydrogens if object contains none."""
    has_h = cmd.count_atoms(f"{obj_name} and elem H") > 0
    if not has_h:
        cmd.h_add(obj_name)

def find_nearest_hydrogen_in_residue(obj_name, atom, max_dist=1.25):
    """
    Find nearest H atom in same OBJECT + chain + resi within max_dist from given atom.
    Does NOT rely on atom.model or atom.segi (which vary across PyMOL builds).
    """
    sele = f"{obj_name} and chain {atom.chain} and resi {atom.resi} and elem H"
    hs = cmd.get_model(sele).atom

    best = None
    best_d = 1e9
    for h in hs:
        d = distance_coords(atom.coord, h.coord)
        if d < max_dist and d < best_d:
            best = h
            best_d = d
    return best


# INTERACTION RULES

def is_hydrophobic(atom1, atom2, dist):
    hydrophobic_elems = {"C", "S"}
    e1 = get_element(atom1)
    e2 = get_element(atom2)
    return dist <= 4.0 and (e1 in hydrophobic_elems) and (e2 in hydrophobic_elems)

def is_hbond(obj1_name, atom1, obj2_name, atom2, dist):
    """
    Symmetric H-bond check:
    - If atom1 has nearby H (in obj1), check angle atom1-H-atom2
    - If atom2 has nearby H (in obj2), check angle atom2-H-atom1
    Distance criterion: dist <= hbond_distance_cutoff
    """
    if dist > hbond_distance_cutoff:
        return False

    da_elems = {"N", "O", "F"}
    e1 = get_element(atom1)
    e2 = get_element(atom2)
    if e1 not in da_elems or e2 not in da_elems:
        return False

    h1 = find_nearest_hydrogen_in_residue(obj1_name, atom1)
    if h1:
        ang = calculate_angle(atom1, h1, atom2)
        if ang >= hbond_angle_cutoff:
            return True

    h2 = find_nearest_hydrogen_in_residue(obj2_name, atom2)
    if h2:
        ang = calculate_angle(atom2, h2, atom1)
        if ang >= hbond_angle_cutoff:
            return True

    return False

def is_salt_bridge(atom1, atom2, dist):
    positive_atoms = {"NZ", "NH1", "NH2", "NE", "ND1", "ND2"}
    negative_atoms = {"OD1", "OD2", "OE1", "OE2"}
    return dist <= 4.0 and (
        (atom1.name in positive_atoms and atom2.name in negative_atoms) or
        (atom1.name in negative_atoms and atom2.name in positive_atoms)
    )

def is_aromatic(atom1, atom2, dist):
    aromatic_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 5.0 and atom1.name in aromatic_atoms and atom2.name in aromatic_atoms

def is_cation_pi(atom1, atom2, dist):
    cation_atoms = {"NZ", "NH1", "NH2", "NE"}
    pi_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 6.0 and (
        (atom1.name in cation_atoms and atom2.name in pi_atoms) or
        (atom2.name in cation_atoms and atom1.name in pi_atoms)
    )

def is_covalent(dist):
    return dist < 1.6

def is_vdw(dist):
    # vdw is assigned only after other rules by ordering in the main loop
    return 1.6 < dist <= 4.0


# MAIN SCRIPT

# Validate objects exist
ensure_object_exists(rec_file)
ensure_object_exists(lig_file)
if rec_file == lig_file:
    raise ValueError("rec and lig must be separate objects (different names).")

# Prepare output path
os.makedirs(output_directory, exist_ok=True)
output_path = os.path.join(output_directory, output_filename)

# Optionally add hydrogens
if ensure_hydrogens:
    ensure_h_atoms(rec_file)
    ensure_h_atoms(lig_file)

interactions = {
    "hydrophobic": [], "hbond": [], "salt_bridge": [], "covalent": [],
    "vdw": [], "aromatic": [], "cation_pi": [], "other": []
}

cmd.feedback("push")
cmd.feedback("disable", "all", "actions")
cmd.feedback("enable", "cmd", "results")

print(f"Running interaction analysis between '{rec_file}' and '{lig_file}'...")

# Select interacting atoms/residues
cmd.select(
    "close_contacts",
    f"({rec_file} within {distance_cutoff} of {lig_file}) or ({lig_file} within {distance_cutoff} of {rec_file})"
)
if cmd.count_atoms("close_contacts") == 0:
    cmd.feedback("pop")
    raise RuntimeError(
        f"No close contacts found within {distance_cutoff} Å.\n"
        "Check that 'rec' and 'lig' are positioned correctly and are separate objects."
    )

cmd.select("interacting_residues", "byres close_contacts")

rec_atoms = cmd.get_model(f"{rec_file} and close_contacts").atom
lig_atoms = cmd.get_model(f"{lig_file} and close_contacts").atom

# Determine interactions
for rec_atom in rec_atoms:
    for lig_atom in lig_atoms:
        dist = cmd.get_distance(
            f"{rec_file}//{rec_atom.chain}/{rec_atom.resi}/{rec_atom.name}",
            f"{lig_file}//{lig_atom.chain}/{lig_atom.resi}/{lig_atom.name}"
        )
        if dist <= distance_cutoff:
            if is_hbond(rec_file, rec_atom, lig_file, lig_atom, dist):
                interactions["hbond"].append((rec_atom, lig_atom, dist))
            elif is_salt_bridge(rec_atom, lig_atom, dist):
                interactions["salt_bridge"].append((rec_atom, lig_atom, dist))
            elif is_hydrophobic(rec_atom, lig_atom, dist):
                interactions["hydrophobic"].append((rec_atom, lig_atom, dist))
            elif is_aromatic(rec_atom, lig_atom, dist):
                interactions["aromatic"].append((rec_atom, lig_atom, dist))
            elif is_cation_pi(rec_atom, lig_atom, dist):
                interactions["cation_pi"].append((rec_atom, lig_atom, dist))
            elif is_covalent(dist):
                interactions["covalent"].append((rec_atom, lig_atom, dist))
            elif is_vdw(dist):
                interactions["vdw"].append((rec_atom, lig_atom, dist))
            else:
                interactions["other"].append((rec_atom, lig_atom, dist))

# Write results to file
with open(output_path, "w", encoding="utf-8") as f:
    for interaction_type, pairs in interactions.items():
        f.write(f"{interaction_type.upper()} INTERACTIONS:\n")
        for atom1, atom2, dist in pairs:
            f.write(
                f"{atom1.chain}/{atom1.resn}`{atom1.resi}/{atom1.name} -- "
                f"{atom2.chain}/{atom2.resn}`{atom2.resi}/{atom2.name} : {dist:.2f} Å\n"
            )
        f.write("\n")

print(f"Results saved to file: {output_path}")

# Visualization
cmd.show("sticks", "interacting_residues")

for interaction_type, color in zip(
    ["hbond", "salt_bridge", "covalent", "hydrophobic"],
    ["cyan", "magenta", "yellow", "orange"]
):
    if interactions[interaction_type]:
        for atom1, atom2, _dist in interactions[interaction_type]:
            cmd.distance(
                f"{interaction_type}_bonds",
                f"{rec_file}//{atom1.chain}/{atom1.resi}/{atom1.name}",
                f"{lig_file}//{atom2.chain}/{atom2.resi}/{atom2.name}"
            )
        cmd.show("dashes", f"{interaction_type}_bonds")
        cmd.color(color, f"{interaction_type}_bonds")

cmd.feedback("pop")
