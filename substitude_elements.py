#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
from aim2dat.strct import StructureCollection
from aim2dat.strct import StructureOperations

# Initialize structure collection and load the Zn_MOF5_H structure
strct = StructureCollection()
strct.append_from_file("Zn_MOF5_H", "Zn_MOF5_H_conv.xyz")

# Initialize StructureOperations with the loaded structure
strct_op = StructureOperations(structures=strct)

# Substitute Zn with Ca and save the structure
strct_subst_ca = strct_op.substitute_elements(
    strct_op.structures.labels,
    [("Zn", "Ca")],
    change_label=True
)

# Verify the label and save if the substitution was successful
if strct_subst_ca.get_structure("Zn_MOF5_H_subst-ZnCa") is not None:
    strct_subst_ca.get_structure("Zn_MOF5_H_subst-ZnCa").to_file("Ca_MOF5_H_conv.xyz")
else:
    print("Substitution for Zn -> Ca was unsuccessful or label is incorrect.")

# Substitute Zn with Ba and save the structure
strct_subst_ba = strct_op.substitute_elements(
    strct_op.structures.labels,
    [("Zn", "Ba")],
    change_label=True
)

# Verify the label and save if the substitution was successful
if strct_subst_ba.get_structure("Zn_MOF5_H_subst-ZnBa") is not None:
    strct_subst_ba.get_structure("Zn_MOF5_H_subst-ZnBa").to_file("Ba_MOF5_H_conv.xyz")
else:
    print("Substitution for Zn -> Ba was unsuccessful or label is incorrect.")

